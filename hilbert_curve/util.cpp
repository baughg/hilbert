#include "util.h"


void write_bitmap(std::string filename, int _width, int _height, int planes,
  uint8_t *dataPtr) {
  FILE *bitmapFile;
  short Type = 19778;
  int Size;
  int Reserved = 0;
  int Offset = 54;
  int headerSize = 40;
  int Width = _width;
  int Height = _height;
  short Planes = 1;
  short BitsPerPixel = 8 * planes;
  int Compression = 0;
  int SizeImage;
  int XPixelsPerMeter = 0;
  int YPixelsPerMeter = 0;
  int ColorsUsed = 0;
  int ColorsImportant = 0;

  int stride = Width * (BitsPerPixel / 8);
  int bytesPerLine = Width * (BitsPerPixel / 8) * sizeof(unsigned char);
  int pad = stride % 4;

  if (pad != 0) {
    stride += (4 - pad);
  }

  SizeImage = stride * Height;
  Size = SizeImage + Offset;
  unsigned char *writeBuffer = new unsigned char[stride];

  bitmapFile = fopen(filename.c_str(), "wb");

  fwrite((const void *)&Type, sizeof(short), 1, bitmapFile);
  fwrite((const void *)&Size, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&Reserved, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&Offset, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&headerSize, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&Width, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&Height, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&Planes, sizeof(short), 1, bitmapFile);
  fwrite((const void *)&BitsPerPixel, sizeof(short), 1, bitmapFile);
  fwrite((const void *)&Compression, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&SizeImage, sizeof(int), 1, bitmapFile);

  fwrite((const void *)&XPixelsPerMeter, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&YPixelsPerMeter, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&ColorsUsed, sizeof(int), 1, bitmapFile);
  fwrite((const void *)&ColorsImportant, sizeof(int), 1, bitmapFile);

  unsigned char *linePtr = dataPtr + (bytesPerLine * (Height - 1));
  //unsigned char *linePtr = dataPtr;

  for (int row = 0; row < Height; row++) {
    memcpy((void *)writeBuffer, (const void *)linePtr, bytesPerLine);
    fwrite((void *)writeBuffer, sizeof(unsigned char), stride, bitmapFile);
    linePtr -= bytesPerLine;
  }

  fclose(bitmapFile);
}