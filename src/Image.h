#pragma once
#include "float3.h"
#include <math.h>

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
    float3 filterAt(float x, float y) const
    {
        const int x0 = (int)floor(x);
        const int y0 = (int)floor(y);
        // we allow ourselves to read past the width/height of the Image because the data is valid
        // and contain the "seamless" data.
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        const float u = float(x - x0);
        const float v = float(y - y0);
        const float one_minus_u = 1 - u;
        const float one_minus_v = 1 - v;
        const float3& c0 = getPixel(x0, y0);
        const float3& c1 = getPixel(x1, y0);
        const float3& c2 = getPixel(x0, y1);
        const float3& c3 = getPixel(x1, y1);
        return c0 * (one_minus_u * one_minus_v) + c1 * (u * one_minus_v) + c2 * (one_minus_u * v) + c3 * (u * v);
    }
    void subset(const Image& image, int x, int y, int w, int h);
};
