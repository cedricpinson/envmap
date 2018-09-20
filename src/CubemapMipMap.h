#pragma once

#include "Cubemap.h"

const size_t CubemapPrefilterMaxLevel = 16;

struct CubemapMipMap
{
    Cubemap levels[CubemapPrefilterMaxLevel];
    int numLevel;
    int padding;
    void init(int num) { numLevel = num; }
};
