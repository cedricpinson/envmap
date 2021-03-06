#pragma once
#include "float3.h"
#include <assert.h>
#include <stdio.h>

// Image is supposed to be 3 rgb float value
struct Image
{
    enum Type { EXR = 0, HDR = 1, RAW = 2 };

    const Image* originalImage = 0;
    int originalX = 0;
    int originalY = 0;

    float3* data;
    int rowInFloat3;
    Type type;
    int width;
    int height;

    float* getData() { return (float*)data; }
    const float* getData() const { return (const float*)data; }
    float3& getPixel(int x, int y) { return data[x + y * rowInFloat3]; }
    const float3& getPixel(int x, int y) const { return data[x + y * rowInFloat3]; }
    inline void filterAt(float* pixel, float x, float y) const
    {
        const float x0 = (int)x; // floor the int
        const float y0 = (int)y;
        // we allow ourselves to read past the width/height of the Image because the data is valid
        // and contain the "seamless" data.
        const float u = x - x0;
        const float v = y - y0;
        const float one_minus_u = 1.f - u;
        const float one_minus_v = 1.f - v;

        const int xInt0 = (int)x0;
        const int yInt0 = (int)y0;

        const float* c0 = data[xInt0 + yInt0 * rowInFloat3]._v; // x0,y0
        const float* c1;                                        // x1,y0

        // it's commented part is equivalent with part that clamp
        // const float3* p1 = &getPixel(xInt0, yInt0 + 1);
        // const float3& c2 = *p1;
        // const float3& c3 = p1[1];
        // const float* c2 = c0 + rowInFloat3 * 3; // x0,y1
        // const float* c3 = c2 + 3;               // x1,y1

        const float* c2; // x0,y1
        const float* c3; // x1,y1
        if (yInt0 == height)
            c2 = c0;
        else
            c2 = c0 + rowInFloat3 * 3; // x0,y1

        if (xInt0 == width)
        {
            c1 = c0; // clamp
            c3 = c2; // clamp
        }
        else
        {
            c1 = c0 + 3;
            c3 = c2 + 3;
        }

        float f1 = (one_minus_u * one_minus_v);
        float f2 = (u * one_minus_v);
        float f3 = (one_minus_u * v);
        float f4 = (u * v);
        pixel[0] = c0[0] * f1 + c1[0] * f2 + c2[0] * f3 + c3[0] * f4;
        pixel[1] = c0[1] * f1 + c1[1] * f2 + c2[1] * f3 + c3[1] * f4;
        pixel[2] = c0[2] * f1 + c1[2] * f2 + c2[2] * f3 + c3[2] * f4;
        // pixel = c0 * (one_minus_u * one_minus_v) + c1 * (u * one_minus_v) + c2 * (one_minus_u * v) + c3 * (u * v);
    }
    void subset(const Image& image, int x, int y, int w, int h);
};
