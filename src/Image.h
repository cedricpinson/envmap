#pragma once
#include "float3.h"

// Image is supposed to be 3 rgb float value
struct Image
{
    enum Type { EXR = 0, HDR = 1, RAW = 2 };

    Type type;
    int rowInFloat3;
    float3* data;
    int width;
    int height;

    float3& getPixel(int x, int y) { return data[x + y * rowInFloat3]; }
    const float3& getPixel(int x, int y) const { return data[x + y * rowInFloat3]; }
    void subset(const Image& image, int x, int y, int w, int h);
};
