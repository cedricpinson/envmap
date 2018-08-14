#include "Cubemap.h"

const size_t CubemapPrefilterMaxLevel = 16;

struct CubemapMipMap
{
    int numLevel;
    Cubemap levels[CubemapPrefilterMaxLevel];
    void init(int num) { numLevel = num; }
};
