#include "HibertCurve.h"



HibertCurve::HibertCurve()
	: order_(0),
	size_(0)
{
}


HibertCurve::~HibertCurve()
{
}


void HibertCurve::curve2d(uint32_t order)
{
	if (!order)
		return;

	order_ = order;

	size_ = 1 << order_;

	uint32_t scale = size_ >> 1;
	uint64_t address;
	uint32_t point_scale = order_ - 1;
	uint32_t point = 0;
	levels_ = 0;

	address_.resize(size_*size_);

	while (scale)
	{
		point = 0;

		for (uint32_t y = 0; y < size_; ++y)
		{
			for (uint32_t x = 0; x < size_; ++x)
			{
				address = 0;
				address = ((y >> point_scale) << 1) | (x >> point_scale);
				address_[point] = (address_[point] << 2) | (address & 0x3ULL);
				point++;
			}
		}

		levels_++;
		scale >>= 1;
	}
}
