#ifndef UTIL_H
#define UTIL_H
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <string>

void write_bitmap(std::string filename, int _width, int _height, int planes,
  uint8_t *dataPtr);

#endif