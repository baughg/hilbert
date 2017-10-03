#ifndef CURVE3D_H
#define CURVE3D_H
#include <stdint.h>

uint32_t mortonToHilbert3D(uint32_t mortonIndex, uint32_t bits);
uint32_t hilbertToMorton3D(uint32_t hilbertIndex, uint32_t bits);
#endif