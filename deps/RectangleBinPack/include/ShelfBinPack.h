/** @file ShelfBinPack.h
	@author Jukka Jylänki

	@brief Implements different bin packer algorithms that use the SHELF data structure.

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include "GuillotineBinPack.h"
#include "Rect.h"

#include <vector>

namespace rbp {

/** ShelfBinPack implements different bin packing algorithms that use the SHELF data structure. ShelfBinPack
also uses GuillotineBinPack for the waste map if it is enabled. */
class ShelfBinPack
{
public:
	/// Default ctor initializes a bin of size (0,0). Call Init() to init an instance.
	ShelfBinPack();

	ShelfBinPack(int width, int height, bool useWasteMap);

	/// Clears all previously packed rectangles and starts packing from scratch into a bin of the given size.
	void Init(int width, int height, bool useWasteMap);

	/// Defines different heuristic rules that can be used in the packing process.
	enum ShelfChoiceHeuristic
	{
		ShelfNextFit, ///< -NF: We always put the new rectangle to the last open shelf.
		ShelfFirstFit, ///< -FF: We test each rectangle against each shelf in turn and pack it to the first where it fits.
		ShelfBestAreaFit, ///< -BAF: Choose the shelf with smallest remaining shelf area.
		ShelfWorstAreaFit, ///< -WAF: Choose the shelf with the largest remaining shelf area.
		ShelfBestHeightFit, ///< -BHF: Choose the smallest shelf (height-wise) where the rectangle fits.
		ShelfBestWidthFit, ///< -BWF: Choose the shelf that has the least remaining horizontal shelf space available after packing.
		ShelfWorstWidthFit, ///< -WWF: Choose the shelf that will have most remainining horizontal shelf space available after packing.
	};

	/// Inserts a single rectangle into the bin. The packer might rotate the rectangle, in which case the returned
	/// struct will have the width and height values swapped.
	/// @param method The heuristic rule to use for choosing a shelf if multiple ones are possible.
	Rect Insert(int width, int height, ShelfChoiceHeuristic method);

	/// Computes the ratio of used surface area to the total bin area.
	float Occupancy() const;

private:
	int binWidth;
	int binHeight;

	/// Stores the starting y-coordinate of the latest (topmost) shelf.
	int currentY;

	/// Tracks the total consumed surface area.
	unsigned long usedSurfaceArea;

	/// If true, the following GuillotineBinPack structure is used to recover the SHELF data structure from losing space.
	bool useWasteMap;
	GuillotineBinPack wasteMap;

	/// Describes a horizontal slab of space where rectangles may be placed.
	struct Shelf
	{
		/// The x-coordinate that specifies where the used shelf space ends.
		/// Space between [0, currentX[ has been filled with rectangles, [currentX, binWidth[ is still available for filling.
		int currentX;

		/// The y-coordinate of where this shelf starts, inclusive.
		int startY;

		/// Specifices the height of this shelf. The topmost shelf is "open" and its height may grow.
		int height;

		/// Lists all the rectangles in this shelf.
		std::vector<Rect> usedRectangles;
	};

	std::vector<Shelf> shelves;

	/// Parses through all rectangles added to the given shelf and adds the gaps between the rectangle tops and the shelf
	/// ceiling into the waste map. This is called only once when the shelf is being closed and a new one is opened.
	void MoveShelfToWasteMap(Shelf &shelf);

	/// Returns true if the rectangle of size width*height fits on the given shelf, possibly rotated.
	/// @param canResize If true, denotes that the shelf height may be increased to fit the object.
	bool FitsOnShelf(const Shelf &shelf, int width, int height, bool canResize) const;

	/// Measures and if desirable, flips width and height so that the rectangle fits the given shelf the best.
	/// @param width [in,out] The width of the rectangle.
	/// @param height [in,out] The height of the rectangle.
	void RotateToShelf(const Shelf &shelf, int &width, int &height) const;

	/// Adds the rectangle of size width*height into the given shelf, possibly rotated.
	/// @param newNode [out] The added rectangle will be returned here.
	void AddToShelf(Shelf &shelf, int width, int height, Rect &newNode);

	/// Returns true if there is still room in the bin to start a new shelf of the given height.
	bool CanStartNewShelf(int height) const;

	/// Creates a new shelf of the given starting height, which will become the topmost 'open' shelf.
	void StartNewShelf(int startingHeight);
};

}
