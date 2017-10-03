#include "HilbertCurve.h"



HilbertCurve::HilbertCurve()
	: order_(0),
	size_(0),
  stride_(0),
  stride_plane_(0)
{
}


HilbertCurve::~HilbertCurve()
{
}


void HilbertCurve::curve2d(uint32_t order)
{
	if (!order)
		return;

	order_ = order;

	size_ = 1 << order_;

	uint32_t scale = size_ >> 1;
	uint64_t address;
  uint64_t address_prev;
	uint32_t point_scale = order_ - 1;
	uint32_t point = 0;
  uint64_t direction;
	levels_ = 0;
  uint32_t x_hat = 0, y_hat = 0;
  uint32_t xa = 0, ya = 0;
	address_.resize(size_*size_);
  direction_.resize(address_.size());

  point = 0;
   
	while (scale)
	{
		point = 0;

		for (uint32_t y = 0; y < size_; ++y)
		{
			for (uint32_t x = 0; x < size_; ++x)
			{
				address = 0;
        xa = x >> point_scale;
        ya = y >> point_scale;
        x_hat = (xa & ~ya) | (xa & ya);
        y_hat = (~xa & ya) | (xa & ~ya);
				address = ((x_hat << 1) & 0x2) | (y_hat & 0x1);
        address &= 0x3ULL;
        address_prev = address_[point] &= 0x3ULL;
        //direction = direction_[point] & 0x3ULL;

        
        address_[point] = (address_[point] << 2) | address;
				point++;
			}
		}
    point_scale--;
		levels_++;
		scale >>= 1;
	}


  point = 0;

  for (uint32_t y = 0; y < size_; ++y)
  {
    for (uint32_t x = 0; x < size_; ++x)
    {
      printf("|%03u,0x%.2X|",point, address_[point]);
      point++;
    }
    printf("\n");
  }
}

void HilbertCurve::swap_state(point_state &a, point_state &b)
{
  point_state a_copy = a;
  a.x = b.x;
  a.y = b.y;
  b.x = a_copy.x;
  b.y = a_copy.y;

  uint32_t point = a.y * stride_ + a.x;
  p_state_[point] = &a;
  point = b.y * stride_ + b.x;
  p_state_[point] = &b;
}

void HilbertCurve::swap_state3d(point_state &a, point_state &b)
{
  point_state a_copy = a;
  a.x = b.x;
  a.y = b.y;
  a.z = b.z;
  b.x = a_copy.x;
  b.y = a_copy.y;
  b.z = a_copy.z;

  uint32_t point = a.z * stride_plane_ + a.y * stride_ + a.x;
  p_state_[point] = &a;
  point = b.z * stride_plane_ + b.y * stride_ + b.x;
  p_state_[point] = &b;
}

void HilbertCurve::grow3d(uint32_t order)
{
  if (!order)
    return;

  order_ = order;
  size_ = 1 << order_;
  stride_ = size_;
  uint32_t offset = 0;
  uint32_t point = 0;
  stride_plane_ = size_ * size_;
  state_.resize(size_*size_*size_);
  p_state_.resize(state_.size());
  uint32_t dim = 0;
  uint32_t max_index = 0;
  uint32_t point_a = 0, point_b = 0;

  for (uint32_t z = 0; z < size_; ++z)
  {
    for (uint32_t y = 0; y < size_; ++y)
    {
      for (uint32_t x = 0; x < size_; ++x)
      {
        state_[point].x = x;
        state_[point].y = y;
        state_[point].z = z;
        p_state_[point] = &state_[point];
        point++;
      }
    }
  }

  state_[0].p_prev = &state_[0]; 
  state_[0].p_next = &state_[stride_plane_];                          // 0 --> 1
  state_[stride_plane_].p_prev = &state_[0];                          // 1 <-- 0
  state_[stride_plane_].p_next = &state_[stride_plane_ + 1];          // 1 --> 2
  state_[stride_plane_ + 1].p_prev = &state_[stride_plane_];          // 2 <-- 1
  state_[stride_plane_ + 1].p_next = &state_[1];                      // 2 --> 3
  state_[1].p_prev = &state_[stride_plane_ + 1];                      // 3 <-- 2
  state_[1].p_next = &state_[stride_ + 1];                            // 3 --> 4
  state_[stride_ + 1].p_prev = &state_[1];                            // 4 <-- 3
  state_[stride_ + 1].p_next = &state_[stride_plane_ + stride_ + 1];  // 4 --> 5
  state_[stride_plane_ + stride_ + 1].p_prev = &state_[stride_ + 1];  // 5 <-- 4
  // 5 --> 6
  state_[stride_plane_ + stride_ + 1].p_next = &state_[stride_plane_ + stride_];
  // 6 <-- 5
  state_[stride_plane_ + stride_].p_prev = &state_[stride_plane_ + stride_ + 1];
  state_[stride_plane_ + stride_].p_next = &state_[stride_];          // 6 --> 7
  state_[stride_].p_prev = &state_[stride_plane_ + stride_];          // 7 <-- 6
  state_[stride_].p_next = NULL;

  // report
  point_state* p_current = &state_[0];
  uint32_t distance = 0;
  uint32_t z_offset = 0;

  int32_t dx = 0;
  int32_t dy = 0;
  int32_t rx = 0;
  int32_t ry = 0;
  point_state* p_qaud_0_start = NULL, *p_qaud_0_end = NULL;
  point_state* p_qaud_1_start = NULL, *p_qaud_1_end = NULL;
  point_state* p_qaud_2_start = NULL, *p_qaud_2_end = NULL;
  point_state* p_qaud_3_start = NULL, *p_qaud_3_end = NULL;

  p_qaud_0_start = &state_[0];

  for (uint32_t o = 2; o <= order_; ++o)
  {
    dim = 1 << (o - 1);
    max_index = dim - 1;
    // copy to quad 1 and 2
    offset = dim * stride_;

    for (uint32_t z = 0; z < dim; ++z)
    {
      z_offset = z * stride_plane_;

      for (uint32_t y = 0; y < dim; ++y)
      {
        for (uint32_t x = y; x < dim; ++x)
        {
          point_a = z_offset + x * stride_ + y;
          point_b = z_offset + y * stride_ + x;

          if (point_b != point_a)
            swap_state3d(*p_state_[point_a], *p_state_[point_b]);
        }
      }
    }
  }

  point = 0;
  while (p_current)
  {
    p_current->distance = distance++;

    printf("%03u. (%02u,%02u,%02u)\n",
      p_current->distance,
      p_current->x,
      p_current->y,
      p_current->z);
    p_current = p_current->p_next;
  }
}

void HilbertCurve::grow2d(uint32_t order)
{
  if (!order)
    return;

  order_ = order;
  size_ = 1 << order_;
  const uint32_t stride = size_;
  uint32_t offset = 0;
  stride_ = stride;
  state_.resize(size_*size_);  
  p_state_.resize(state_.size());
  uint32_t dim = 0;
  uint32_t max_index = 0;
  uint32_t point = 0;
  uint32_t point_b = 0;
  uint32_t point_c = 0;
  uint32_t point_b2 = 0;
  uint32_t point_c2 = 0;
  uint32_t point_b3 = 0;
  uint32_t point_c3 = 0;

  for (uint32_t y = 0; y < size_; ++y)
  {
    for (uint32_t x = 0; x < size_; ++x)
    {
      state_[point].x = x;
      state_[point].y = y; 
      p_state_[point] = &state_[point];
      point++;
    }
  }

  state_[0].p_prev = &state_[0];
  state_[0].p_next = &state_[stride];
  state_[stride].p_prev = &state_[0];
  state_[stride].p_next = &state_[stride+1];
  state_[stride + 1].p_prev = &state_[stride];
  state_[stride + 1].p_next = &state_[1];
  state_[1].p_prev = &state_[stride + 1];
  int32_t dx = 0;
  int32_t dy = 0;
  int32_t rx = 0;
  int32_t ry = 0;
  point_state* p_qaud_0_start = NULL, *p_qaud_0_end = NULL;
  point_state* p_qaud_1_start = NULL, *p_qaud_1_end = NULL;
  point_state* p_qaud_2_start = NULL, *p_qaud_2_end = NULL;
  point_state* p_qaud_3_start = NULL, *p_qaud_3_end = NULL;
  
  p_qaud_0_start = &state_[0];

  for (uint32_t o = 2; o <= order_; ++o)
  {
    dim = 1 << (o-1);
    max_index = dim - 1;
    // copy to quad 1 and 2
    offset = dim * stride;
    
    for (uint32_t y = 0; y < dim; ++y)
    {
      for (uint32_t x = 0; x < dim; ++x)
      {
        point = y * stride + x;
        point_b = point + offset;
        point_b2 = point_b + dim;
        rx = max_index - y;
        ry = max_index - x;
        rx += dim;
        point_b3 = ry * stride + rx;
        

        p_state_[point_b]->x = p_state_[point]->x;
        p_state_[point_b]->y = p_state_[point]->y + dim;

        p_state_[point_b2]->x = p_state_[point_b]->x + dim;
        p_state_[point_b2]->y = p_state_[point_b]->y;

        p_state_[point_b3]->x = rx;
        p_state_[point_b3]->y = ry;

        if (!p_state_[point]->p_next) {
          p_qaud_0_end = p_state_[point];
          p_qaud_1_end = p_state_[point_b];
          p_qaud_2_end = p_state_[point_b2];
          p_qaud_3_end = p_state_[point_b3];
          continue;
        }

        if (p_state_[point]->p_prev == p_state_[point])
        {
          p_qaud_1_start = p_state_[point_b];
          p_qaud_2_start = p_state_[point_b2];
          p_qaud_3_start = p_state_[point_b3];
        }

        // quad 1
        dx = p_state_[point]->p_next->x - p_state_[point]->x;
        dy = p_state_[point]->p_next->y - p_state_[point]->y;
        point_c = (p_state_[point_b]->y + dy) * stride + (p_state_[point_b]->x + dx);
        p_state_[point_b]->p_next = p_state_[point_c];
        p_state_[point_c]->p_prev = p_state_[point_b];
        // quad 2        
        point_c2 = (p_state_[point_b2]->y + dy) * stride + (p_state_[point_b2]->x + dx);
        p_state_[point_b2]->p_next = p_state_[point_c2];
        p_state_[point_c2]->p_prev = p_state_[point_b2];
        // quad 3 
        rx = (max_index - p_state_[point]->p_next->y);
        ry = (max_index - p_state_[point]->p_next->x);
        rx += dim;
        point_c3 = ry * stride + rx;
        p_state_[point_b3]->p_next = p_state_[point_c3];
        p_state_[point_c3]->p_prev = p_state_[point_b3];
      }
    }
    // connect quads
    p_qaud_1_start->p_prev = p_qaud_0_end; 
    p_qaud_0_end->p_next = p_qaud_1_start;

    p_qaud_2_start->p_prev = p_qaud_1_end;
    p_qaud_1_end->p_next = p_qaud_2_start;

    p_qaud_3_start->p_prev = p_qaud_2_end;
    p_qaud_2_end->p_next = p_qaud_3_start;

    // quad 0
    for (uint32_t y = 0; y < dim; ++y)
    {
      for (uint32_t x = y; x < dim; ++x)
      {
        point = x * stride + y;
        point_b = y * stride + x;
        if (point_b != point)
          swap_state(*p_state_[point], *p_state_[point_b]);
      }
    }
  }

  point_state* p_current = p_qaud_0_start;
  uint32_t distance = 0;
  point2d_.resize(state_.size());
  point1d_.resize(state_.size());
  point = 0;
  while (p_current)
  {
    p_current->distance = distance++;
    point2d_[point].x = p_current->x;
    point2d_[point].y = p_current->y;    
    point1d_[stride_*point2d_[point].y + point2d_[point].x] = p_current->distance;
    point++;
    /*printf("%03u. (%02u,%02u)\n", 
      p_current->distance,
      p_current->x,
      p_current->y);*/
    p_current = p_current->p_next;
  }
}

void HilbertCurve::get_xy(const uint32_t &p, uint32_t &x, uint32_t &y)
{
  x = point2d_[p].x;
  y = point2d_[p].y;
}