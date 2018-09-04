#pragma once

#include "Image.h"
#include "float3.h"
#include <math.h>

struct Cubemap
{

    Image faces[6];
    Image image;
    int size;
    int padding;

    // order in opengl
    enum Face {
        PX = 0, // right           +----+
        NX,     // left            | PY |
        PY,     // top        +----+----+----+----+
        NY,     // bottom     | NX | PZ | PX | NZ |
        PZ,     // front      +----+----+----+----+
        NZ      // back            | NY |
                //                 +----+
    };

    struct Address
    {
        Face face;
        float s = 0;
        float t = 0;
    };


    void setImageForFace(Face face, const Image& image);
    void getDirectionFor(float* direction, Face face, int x, int y) const;
    void getDirectionFor(float* direction, Face face, float x, float y) const;
    void setAllFacesFromCross(const Image& image);

    void makeSeamless();
};

inline void Cubemap::getDirectionFor(float* direction, Face face, int x, int y) const
{
    getDirectionFor(direction, face, x + 0.5f, y + 0.5f);
}

inline void Cubemap::getDirectionFor(float* direction, Face face, float x, float y) const
{
    const float scale = 2.0f / size;
    // map [0, dim] to [-1,1] with (-1,-1) at bottom left
    float cx = (x * scale) - 1;
    float cy = 1 - (y * scale);

    const float l = sqrt(cx * cx + cy * cy + 1);
    switch (face)
    {
    case PX:
        direction[0] = 1;
        direction[1] = cy;
        direction[2] = -cx;
        break;
    case NX:
        direction[0] = -1;
        direction[1] = cy;
        direction[2] = cx;
        break;
    case PY:
        direction[0] = cx;
        direction[1] = 1;
        direction[2] = -cy;
        break;
    case NY:
        direction[0] = cx;
        direction[1] = -1;
        direction[2] = cy;
        break;
    case PZ:
        direction[0] = cx;
        direction[1] = cy;
        direction[2] = 1;
        break;
    case NZ:
        direction[0] = -cx;
        direction[1] = cy;
        direction[2] = -1;
        break;
    }

    direction[0] *= 1 / l;
    direction[1] *= 1 / l;
    direction[2] *= 1 / l;
}
