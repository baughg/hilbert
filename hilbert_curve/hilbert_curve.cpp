// hilbert_curve.cpp : Defines the entry point for the console application.
//
#include "util.h"
#include "curve3d.h"
#include "HilbertCurve.h"
#include <stdlib.h>
#include <iostream>
#include <vector>

void d2xy(int n, int d, int &x, int &y);
int i4_power(int i, int j);
void rot(int n, int &x, int &y, int rx, int ry);
int xy2d(int n, int x, int y);

using namespace std;


int main()
{
  const int m = 5; // the index of the Hilbert curve.
  uint32_t h = 0;
  HilbertCurve curve;
  curve.grow2d(m);
  //return 0;

   
  const int d = 0;
  uint32_t x = 0, y = 0;
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
  std::vector<uint32_t> colour_lookup;
  colour_lookup.resize(1 << 24);
  int x_inc = 1;
  int y_inc = 1;
  uint8_t* p_colour = (uint8_t*)&colour_lookup[0];
  int xc = 0;
  int yc = 0;

  for (int z = 0; z < (1 << 8); ++z)
  {
    for (int p = 0; p < (1 << 16); ++p)
    {
      // work
      *p_colour++ = (uint8_t)(xc & 0xff);
      *p_colour++ = (uint8_t)(yc & 0xff);
      *p_colour++ = (uint8_t)(z & 0xff);
      p_colour++;
      xc += x_inc;
      if (xc >= 256)
      {
        yc += y_inc;
        xc = 255;
        x_inc = -1;
      }
      else if (xc <= -1)
      {
        yc += y_inc;
        xc = 0;
        x_inc = 1;
      }
    }

    if (yc >= 256)
    {     
      yc = 255;
      y_inc = -1;
    }
    else if (yc <= -1)
    {      
      yc = 0;
      y_inc = 1;
    }
  }
#define d2xy_run
#ifdef d2xy_run
  memset(p_first_pixel, 0xff, points * 3);

  for (int p = 0; p < points; ++p) {
    curve.get_xy(p, x, y);
    //d2xy(m, p, x, y);
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