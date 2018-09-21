#pragma once

#include <stdint.h>

void encodeLUV(uint8_t* luvDst, const float* rgbSrc);
void decodeLUV(float* rgb, const uint8_t* luv);

void encodeRGBM(uint8_t* rgbmDst, const float* rgbSrc);
void decodeRGBM(float* rgb, const uint8_t* rgbm);


void tonemap(uint8_t* rgbDst, const float* rgbSrc);
