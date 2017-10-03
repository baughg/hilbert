#ifndef HILBERT_CURVE_H
#define HILBERT_CURVE_H
#include <stdint.h>
#include <vector>

typedef struct point_state
{
  point_state()
  {    
    distance = 0;    
    p_next = 0;
    p_prev = 0;
    x = 0;
    y = 0;
    z = 0;
  }
  uint32_t x;
  uint32_t y;
  uint32_t z;
  uint32_t distance;  
  point_state* p_next;
  point_state* p_prev;
}point_state;

typedef struct point2d
{
  point2d()
  {
    x = 0;
    y = 0;
  }
  uint32_t x;
  uint32_t y;
}point2d;

typedef struct point3d
{
  point3d()
  {
    x = 0;
    y = 0;
    z = 0;
    order = 0;
  }

  uint32_t x;
  uint32_t y;
  uint32_t z;
  uint32_t order;
}point3d;

class HilbertCurve
{
public:
  HilbertCurve();
	~HilbertCurve();
	void curve2d(uint32_t order);
  void grow2d(uint32_t order);
  void grow3d(uint32_t order);
  void get_xy(const uint32_t &p, uint32_t &x, uint32_t &y);
private:
  void swap_state(point_state &a, point_state &b);
  void swap_state3d(point_state &a, point_state &b);
	uint32_t order_;
	std::vector<uint64_t> address_;
  std::vector<uint64_t> direction_;
  std::vector<point_state> state_;
  std::vector<point_state*> p_state_;
  std::vector<point2d> point2d_;
  std::vector<uint32_t> point1d_;
	uint32_t point_count_;
	uint32_t size_;
	uint32_t levels_;
  uint32_t stride_;
  uint32_t stride_plane_;
};
#endif

