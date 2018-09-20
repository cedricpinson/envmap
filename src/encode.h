#pragma once

#include <stdint.h>

void encodeLUV(uint8_t* luvDst, const float* rgbSrc);
void decodeLUV(float* rgb, const uint8_t* luv);

void encodeRGBM(float rgb[3], uint8_t rgbm[4]);
void decodeRGBM(uint8_t rgbm[4], float rgb[3]);

void tonemap(uint8_t* rgbDst, const float* rgbSrc);
