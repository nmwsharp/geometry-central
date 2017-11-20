/** @file ShelfBinPack.cpp
	@author Jukka Jylänki

	@brief Implements different bin packer algorithms that use the SHELF data structure.

	This work is released to Public Domain, do whatever you want with it.
*/
#include <utility>
#include <iostream>

#include <cassert>
#include <cstring>

#include "ShelfBinPack.h"

namespace rbp {

using namespace std;

ShelfBinPack::ShelfBinPack()
:binWidth(0),
binHeight(0),
currentY(0),
usedSurfaceArea(0),
useWasteMap(false)
{
}

ShelfBinPack::ShelfBinPack(int width, int height, bool useWasteMap)
{
	Init(width, height, useWasteMap);
}

void ShelfBinPack::Init(int width, int height, bool useWasteMap_)
{
	useWasteMap = useWasteMap_;
	binWidth = width;
	binHeight = height;

	currentY = 0;
	usedSurfaceArea = 0;

	shelves.clear();
	StartNewShelf(0);

	if (useWasteMap)
	{
		wasteMap.Init(width, height);
		wasteMap.GetFreeRectangles().clear();
	}
}

bool ShelfBinPack::CanStartNewShelf(int height) const
{
	return shelves.back().startY + shelves.back().height + height <= binHeight;
}

void ShelfBinPack::StartNewShelf(int startingHeight)
{
	if (shelves.size() > 0)
	{
		assert(shelves.back().height != 0);
		currentY += shelves.back().height;
		
		assert(currentY < binHeight);
	}

	Shelf shelf;
	shelf.currentX = 0;
	shelf.height = startingHeight;
	shelf.startY = currentY;

	assert(shelf.startY + shelf.height <= binHeight);
	shelves.push_back(shelf);
}

bool ShelfBinPack::FitsOnShelf(const Shelf &shelf, int width, int height, bool canResize) const
{
	const int shelfHeight = canResize ? (binHeight - shelf.startY) : shelf.height;
	if ((shelf.currentX + width <= binWidth && height <= shelfHeight) ||
		(shelf.currentX + height <= binWidth && width <= shelfHeight))
		return true;
	else
		return false;
}

void ShelfBinPack::RotateToShelf(const Shelf &shelf, int &width, int &height) const
{	
	// If the width > height and the long edge of the new rectangle fits vertically onto the current shelf,
	// flip it. If the short edge is larger than the current shelf height, store
	// the short edge vertically.
	if ((width > height && width > binWidth - shelf.currentX) ||
		(width > height && width < shelf.height) ||
		(width < height && height > shelf.height && height <= binWidth - shelf.currentX))
		swap(width, height);
}

void ShelfBinPack::AddToShelf(Shelf &shelf, int width, int height, Rect &newNode)
{
	assert(FitsOnShelf(shelf, width, height, true));

	// Swap width and height if the rect fits better that way.
	RotateToShelf(shelf, width, height);

	// Add the rectangle to the shelf.
	newNode.x = shelf.currentX;
	newNode.y = shelf.startY;
	newNode.width = width;
	newNode.height = height;
	shelf.usedRectangles.push_back(newNode);

	// Advance the shelf end position horizontally.
	shelf.currentX += width;
	assert(shelf.currentX <= binWidth);

	// Grow the shelf height.
	shelf.height = max(shelf.height, height);
	assert(shelf.height <= binHeight);

	usedSurfaceArea += width * height;
}

Rect ShelfBinPack::Insert(int width, int height, ShelfChoiceHeuristic method)
{
	Rect newNode;

	// First try to pack this rectangle into the waste map, if it fits.
	if (useWasteMap)
	{
		newNode = wasteMap.Insert(width, height, true, GuillotineBinPack::RectBestShortSideFit, 
			GuillotineBinPack::SplitMaximizeArea);
		if (newNode.height != 0)
		{
			// Track the space we just used.
			usedSurfaceArea += width * height;

			return newNode;
		}
	}

	switch(method)
	{
	case ShelfNextFit:		
		if (FitsOnShelf(shelves.back(), width, height, true))
		{
			AddToShelf(shelves.back(), width, height, newNode);
			return newNode;
		}
		break;
	case ShelfFirstFit:		
		for(size_t i = 0; i < shelves.size(); ++i)
			if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
			{
				AddToShelf(shelves[i], width, height, newNode);
				return newNode;
			}
		break;

	case ShelfBestAreaFit:
		{
			// Best Area Fit rule: Choose the shelf with smallest remaining shelf area.
			Shelf *bestShelf = 0;
			unsigned long bestShelfSurfaceArea = (unsigned long)-1;
			for(size_t i = 0; i < shelves.size(); ++i)
			{
				// Pre-rotate the rect onto the shelf here already so that the area fit computation
				// is done correctly.
				RotateToShelf(shelves[i], width, height);
				if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
				{
					unsigned long surfaceArea = (binWidth - shelves[i].currentX) * shelves[i].height;
					if (surfaceArea < bestShelfSurfaceArea)
					{
						bestShelf = &shelves[i];
						bestShelfSurfaceArea = surfaceArea;
					}
				}
			}

			if (bestShelf)
			{
				AddToShelf(*bestShelf, width, height, newNode);
				return newNode;
			}
		}
		break;

	case ShelfWorstAreaFit:
		{
			// Worst Area Fit rule: Choose the shelf with smallest remaining shelf area.
			Shelf *bestShelf = 0;
			int bestShelfSurfaceArea = -1;
			for(size_t i = 0; i < shelves.size(); ++i)
			{
				// Pre-rotate the rect onto the shelf here already so that the area fit computation
				// is done correctly.
				RotateToShelf(shelves[i], width, height);
				if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
				{
					int surfaceArea = (binWidth - shelves[i].currentX) * shelves[i].height;
					if (surfaceArea > bestShelfSurfaceArea)
					{
						bestShelf = &shelves[i];
						bestShelfSurfaceArea = surfaceArea;
					}
				}
			}

			if (bestShelf)
			{
				AddToShelf(*bestShelf, width, height, newNode);
				return newNode;
			}
		}
		break;

	case ShelfBestHeightFit:
		{
			// Best Height Fit rule: Choose the shelf with best-matching height.
			Shelf *bestShelf = 0;
			int bestShelfHeightDifference = 0x7FFFFFFF;
			for(size_t i = 0; i < shelves.size(); ++i)
			{
				// Pre-rotate the rect onto the shelf here already so that the height fit computation
				// is done correctly.
				RotateToShelf(shelves[i], width, height);
				if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
				{
					int heightDifference = max(shelves[i].height - height, 0);
					assert(heightDifference >= 0);

					if (heightDifference < bestShelfHeightDifference)
					{
						bestShelf = &shelves[i];
						bestShelfHeightDifference = heightDifference;
					}
				}
			}

			if (bestShelf)
			{
				AddToShelf(*bestShelf, width, height, newNode);
				return newNode;
			}
		}
		break;

	case ShelfBestWidthFit:
		{
			// Best Width Fit rule: Choose the shelf with smallest remaining shelf width.
			Shelf *bestShelf = 0;
			int bestShelfWidthDifference = 0x7FFFFFFF;
			for(size_t i = 0; i < shelves.size(); ++i)
			{
				// Pre-rotate the rect onto the shelf here already so that the height fit computation
				// is done correctly.
				RotateToShelf(shelves[i], width, height);
				if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
				{
					int widthDifference = binWidth - shelves[i].currentX - width;
					assert(widthDifference >= 0);

					if (widthDifference < bestShelfWidthDifference)
					{
						bestShelf = &shelves[i];
						bestShelfWidthDifference = widthDifference;
					}
				}
			}

			if (bestShelf)
			{
				AddToShelf(*bestShelf, width, height, newNode);
				return newNode;
			}
		}
		break;

	case ShelfWorstWidthFit:
		{
			// Worst Width Fit rule: Choose the shelf with smallest remaining shelf width.
			Shelf *bestShelf = 0;
			int bestShelfWidthDifference = -1;
			for(size_t i = 0; i < shelves.size(); ++i)
			{
				// Pre-rotate the rect onto the shelf here already so that the height fit computation
				// is done correctly.
				RotateToShelf(shelves[i], width, height);
				if (FitsOnShelf(shelves[i], width, height, i == shelves.size()-1))
				{
					int widthDifference = binWidth - shelves[i].currentX - width;
					assert(widthDifference >= 0);

					if (widthDifference > bestShelfWidthDifference)
					{
						bestShelf = &shelves[i];
						bestShelfWidthDifference = widthDifference;
					}
				}
			}

			if (bestShelf)
			{
				AddToShelf(*bestShelf, width, height, newNode);
				return newNode;
			}
		}
		break;

	}

	// The rectangle did not fit on any of the shelves. Open a new shelf.

	// Flip the rectangle so that the long side is horizontal.
	if (width < height && height <= binWidth)
		swap(width, height);

	if (CanStartNewShelf(height))
	{
		if (useWasteMap)
			MoveShelfToWasteMap(shelves.back());
		StartNewShelf(height);
		assert(FitsOnShelf(shelves.back(), width, height, true));
		AddToShelf(shelves.back(), width, height, newNode);
		return newNode;
	}
/*
	///\todo This is problematic: If we couldn't start a new shelf - should we give up
	///      and move all the remaining space of the bin for the waste map to track,
	///      or should we just wait if the next rectangle would fit better? For now,
	///      don't add the leftover space to the waste map. 
	else if (useWasteMap)
	{
		assert(binHeight - shelves.back().startY >= shelves.back().height);
		shelves.back().height = binHeight - shelves.back().startY;
		if (shelves.back().height > 0)
			MoveShelfToWasteMap(shelves.back());

		// Try to pack the rectangle again to the waste map.
		GuillotineBinPack::Node node = wasteMap.Insert(width, height, true, 1, 3);
		if (node.height != 0)
		{
			newNode.x = node.x;
			newNode.y = node.y;
			newNode.width = node.width;
			newNode.height = node.height;
			return newNode;
		}
	}
*/

	// The rectangle didn't fit.
	memset(&newNode, 0, sizeof(Rect));
	return newNode;
}

void ShelfBinPack::MoveShelfToWasteMap(Shelf &shelf)
{
	std::vector<Rect> &freeRects = wasteMap.GetFreeRectangles();

	// Add the gaps between each rect top and shelf ceiling to the waste map.
	for(size_t i = 0; i < shelf.usedRectangles.size(); ++i)
	{
		const Rect &r = shelf.usedRectangles[i];
		Rect newNode;
		newNode.x = r.x;
		newNode.y = r.y + r.height;
		newNode.width = r.width;
		newNode.height = shelf.height - r.height;
		if (newNode.height > 0)
			freeRects.push_back(newNode);		
	}
	shelf.usedRectangles.clear();

	// Add the space after the shelf end (right side of the last rect) and the shelf right side. 
	Rect newNode;
	newNode.x = shelf.currentX;
	newNode.y = shelf.startY;
	newNode.width = binWidth - shelf.currentX;
	newNode.height = shelf.height;
	if (newNode.width > 0)
		freeRects.push_back(newNode);

	// This shelf is DONE.
	shelf.currentX = binWidth;

	// Perform a rectangle merge step.
	wasteMap.MergeFreeList();
}

/// Computes the ratio of used surface area to the bin area.
float ShelfBinPack::Occupancy() const
{
	return (float)usedSurfaceArea / (binWidth * binHeight);
}

}