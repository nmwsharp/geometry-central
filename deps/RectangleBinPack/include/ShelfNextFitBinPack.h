/** @file ShelfNextFitBinPack.h
	@author Jukka Jylänki

	@brief Implements the naive Shelf Next Fit bin packer algorithm.

	This algorithm is not recommended for real use at all - its only advantage is that it
	consumes only a constant amount of memory, whereas the other packers in this library
	use at least a linear amount of memory.

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include <vector>

namespace rbp {

class ShelfNextFitBinPack
{
public:
	struct Node
	{
		int x;
		int y;
		int width;
		int height;

		bool flipped;
	};

	void Init(int width, int height);

	Node Insert(int width, int height);

	/// Computes the ratio of used surface area.
	float Occupancy() const;

private:
	int binWidth;
	int binHeight;

	int currentX;
	int currentY;
	int shelfHeight;

	unsigned long usedSurfaceArea;
};

}
