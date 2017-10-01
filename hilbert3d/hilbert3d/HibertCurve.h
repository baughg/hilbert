#ifndef HILBERT_CURVE_H
#define HILBERT_CURVE_H
#include <stdint.h>
#include <vector>

class HibertCurve
{
public:
	HibertCurve();
	~HibertCurve();
	void curve2d(uint32_t order);
private:
	uint32_t order_;
	std::vector<uint64_t> address_;
	uint32_t point_count_;
	uint32_t size_;
	uint32_t levels_;
};
#endif

