// hilbert_curve.cpp : Defines the entry point for the console application.
//
#include "util.h"
#include "curve3d.h"
#include "HilbertCurve.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>

void d2xy(int n, int d, int &x, int &y);
int i4_power(int i, int j);
void rot(int n, int &x, int &y, int rx, int ry);
int xy2d(int n, int x, int y);

using namespace std;

bool sort_hilbert(const point3d &lhs, const point3d &rhs)
{
  return lhs.order < rhs.order;
}

uint32_t SeparateBy1(unsigned int x) {
  x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x << 8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x << 4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x << 2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x << 1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  return x;
}


uint32_t MortonCode2(unsigned int x, unsigned int y) {
  return SeparateBy1(x) | (SeparateBy1(y) << 1);
}

uint32_t MortonCode3(unsigned int x, unsigned int y, unsigned int z) {
  uint32_t out = 0;
  uint32_t x_lsb, y_lsb, z_lsb, lsb;

  for (int b = 0; b < 8; ++b)
  {
    x_lsb = x & 0x1;
    y_lsb = y & 0x1;
    z_lsb = z & 0x1;
    x >>= 1;
    y >>= 1;
    z >>= 1;

    lsb = x_lsb;
    lsb |= (y_lsb << 1);
    lsb |= (z_lsb << 2);
    lsb <<= (b * 3);
    out |= lsb;
  }

  return out;
}

int main()
{
  const int m = 12; // the index of the Hilbert curve.
	const int rgb_bits = 8;
  uint32_t h = 0;
  //HilbertCurve curve;

  std::vector<uint32_t> colour_lookup;
  colour_lookup.resize(1 << 24);
  int x_inc = 1;
  int y_inc = 1;
  uint8_t* p_colour = (uint8_t*)&colour_lookup[0];

  FILE* colour_file = NULL;

#define GENERATE_COLOUR_LUT
#ifdef GENERATE_COLOUR_LUT
  //curve.grow2d(m);
  //curve.grow3d(2);
  uint32_t limit = (1 << rgb_bits);
  uint32_t hcode = 0;
  uint32_t code = 0;
  std::vector<point3d> point3d_ordered;
  uint32_t colour_count_all = limit * limit * limit;
  //colour_count_all *= (limit*limit);
  point3d_ordered.resize(colour_count_all);

  for(int z = 0; z < limit; ++z)
    for (int y = 0; y < limit; ++y)
      for (int x = 0; x < limit; ++x) {        
        code = x + (y * limit) + (z * limit * limit);
        code = MortonCode3(x, y, z);
        hcode = mortonToHilbert3D(code, rgb_bits);
        point3d_ordered[code].order = hcode;
        point3d_ordered[code].x = x;
        point3d_ordered[code].y = y;
        point3d_ordered[code].z = z;
        //printf("%02u. (%02u,%02u,%02u)\n", hcode,x,y,z);
      }

  std::sort(point3d_ordered.begin(), point3d_ordered.end(), sort_hilbert);
  size_t point_orders = point3d_ordered.size();
  
  fopen_s(&colour_file, "rgb_line.dat", "wb");
  for (size_t p = 0; p < point_orders; ++p)
  {
    colour_lookup[p] = point3d_ordered[p].x + (point3d_ordered[p].y * 256);
    colour_lookup[p] += (point3d_ordered[p].z * 256 * 256);
    /*printf("%02u. (%02u,%02u,%02u)\n", 
      point3d_ordered[p].order, 
      point3d_ordered[p].x, 
      point3d_ordered[p].y, 
      point3d_ordered[p].z);*/
  }
  
  if (colour_file)
  {
    fwrite(&colour_lookup[0], sizeof(uint32_t), colour_lookup.size(), colour_file);
    fclose(colour_file);
  }
#else
  fopen_s(&colour_file, "rgb_line.dat", "rb");

  if (colour_file)
  {
    fread(&colour_lookup[0], sizeof(uint32_t), colour_lookup.size(), colour_file);
    fclose(colour_file);
  }
  else
    exit(0);
#endif
   
  const int d = 0;
  int x = 0, y = 0;
  int n = 1 << m;
  int points = n * n;
  uint32_t colour_step = 1 << 24;
  colour_step /= points;
  std::vector<uint8_t> hilbert_image;
  hilbert_image.resize(points * 3);
  uint8_t* p_pixel = (uint8_t*)&hilbert_image[0];
  uint8_t* p_first_pixel = p_pixel;
  uint32_t colour = 0;
  uint32_t channel = 0;
  uint32_t channel2 = 0;
  char file_out_name[256];
  
  //int xc = 0;
  //int yc = 0;

  //for (int z = 0; z < (1 << 8); ++z)
  //{
  //  for (int p = 0; p < (1 << 16); ++p)
  //  {
  //    // work
  //    *p_colour++ = (uint8_t)(xc & 0xff);
  //    *p_colour++ = (uint8_t)(yc & 0xff);
  //    *p_colour++ = (uint8_t)(z & 0xff);
  //    p_colour++;
  //    xc += x_inc;
  //    if (xc >= 256)
  //    {
  //      yc += y_inc;
  //      xc = 255;
  //      x_inc = -1;
  //    }
  //    else if (xc <= -1)
  //    {
  //      yc += y_inc;
  //      xc = 0;
  //      x_inc = 1;
  //    }
  //  }

  //  if (yc >= 256)
  //  {     
  //    yc = 255;
  //    y_inc = -1;
  //  }
  //  else if (yc <= -1)
  //  {      
  //    yc = 0;
  //    y_inc = 1;
  //  }
  //}
//#define d2xy_run
#ifdef d2xy_run
  memset(p_first_pixel, 0xff, points * 3);

  for (int p = 0; p < points; ++p) {
    //curve.get_xy(p, x, y);
    d2xy(m, p, x, y);
    //printf("%03d. (%03d,%03d)\n",p,x,y);
    //channel = colour;
    /*channel = p % 256;
    channel += 128;
    channel |= (channel << 8);
    channel |= (channel << 8);*/
    channel = 0;
    p_pixel = p_first_pixel;
    p_pixel += ((y*n * 3) + (x * 3));

    p_pixel[2] = (uint8_t)(channel & 0xff);
    channel >>= 8;
    p_pixel[1] = (uint8_t)(channel & 0xff);
    channel >>= 8;
    p_pixel[0] = (uint8_t)(channel & 0xff); 
    
    sprintf(file_out_name,"trace\\point%05d.bmp",p);
    write_bitmap(file_out_name, n, n, 3, p_first_pixel);
    memset(p_first_pixel, 0xff, points * 3);

    //colour += colour_step;
  }
#else
  for (y = 0; y < n; ++y)
  {
    for (x = 0; x < n; ++x)
    {
      colour = xy2d(m, x, y);
      /*channel = colour % 256;*/
      channel = colour_lookup[colour];
      //channel = mortonToHilbert3D(colour, 24);
      //channel %= (1 << 24);
      //channel %= 256;
      //channel += 128;
      //channel |= (channel << 8);
      //channel |= (channel << 8);
      //channel &= 0x00ffff00;
      //channel |= (channel2 & 0xff);
      p_pixel[2] = (uint8_t)(channel & 0xff);
      channel >>= 8;
      p_pixel[1] = (uint8_t)(channel & 0xff);
      channel >>= 8;
      p_pixel[0] = (uint8_t)(channel & 0xff);
      p_pixel += 3;
      //colour += colour_step;
    }
  }
  write_bitmap("test2.bmp", n, n, 3, &hilbert_image[0]);
#endif
  
  return 0;
}

//****************************************************************************80

void d2xy(int m, int d, int &x, int &y)

//****************************************************************************80
//
//  Purpose:
//
//    D2XY converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, int D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
//    Output, int &X, &Y, the Cartesian coordinates of the cell.
//    0 <= X, Y < N.
//
{
  int n;
  int rx;
  int ry;
  int s;
  int t = d;

  n = i4_power(2, m);

  x = 0;
  y = 0;
  for (s = 1; s < n; s = s * 2)
  {
    rx = (1 & (t / 2));
    ry = (1 & (t ^ rx));
    rot(s, x, y, rx, ry);
    x = x + s * rx;
    y = y + s * ry;
    t = t / 4;
  }
  return;
}
//****************************************************************************80

int i4_power(int i, int j)

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if (j < 0)
  {
    if (i == 1)
    {
      value = 1;
    }
    else if (i == 0)
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit(1);
    }
    else
    {
      value = 0;
    }
  }
  else if (j == 0)
  {
    if (i == 0)
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit(1);
    }
    else
    {
      value = 1;
    }
  }
  else if (j == 1)
  {
    value = i;
  }
  else
  {
    value = 1;
    for (k = 1; k <= j; k++)
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void rot(int n, int &x, int &y, int rx, int ry)

//****************************************************************************80
//
//  Purpose:
//
//    ROT rotates and flips a quadrant appropriately
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int N, the length of a side of the square.  N must be a power of 2.
//
//    Input/output, int &X, &Y, the input and output coordinates of a point.
//
//    Input, int RX, RY, ???
//
{
  int t;

  if (ry == 0)
  {
    //
    //  Reflect.
    //
    if (rx == 1)
    {
      x = n - 1 - x;
      y = n - 1 - y;
    }
    //
    //  Flip.
    //
    t = x;
    x = y;
    y = t;
  }
  return;
}
//****************************************************************************80


//****************************************************************************80

int xy2d(int m, int x, int y)

//****************************************************************************80
//
//  Purpose:
//
//    XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
//
//  Discussion:
//
//    It is assumed that a square has been divided into an NxN array of cells,
//    where N is a power of 2.
//
//    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
//    right corner.
//
//  Modified:
//
//    24 December 2015
//
//  Parameters:
//
//    Input, int M, the index of the Hilbert curve.
//    The number of cells is N=2^M.
//    0 < M.
//
//    Input, int X, Y, the Cartesian coordinates of a cell.
//    0 <= X, Y < N.
//
//    Output, int XY2D, the Hilbert coordinate of the cell.
//    0 <= D < N * N.
//
{
  int d = 0;
  int n;
  int rx;
  int ry;
  int s;

  n = i4_power(2, m);

  for (s = n / 2; s > 0; s = s / 2)
  {
    rx = (x & s) > 0;
    ry = (y & s) > 0;
    d = d + s * s * ((3 * rx) ^ ry);
    rot(s, x, y, rx, ry);
  }
  return d;
}