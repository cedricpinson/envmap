#pragma once

#include "Image.h"
#include "float3.h"
#include <math.h>

struct Cubemap
{

    Image image;
    int size;

    enum Face {
        NX = 0, // left            +----+
        PX,     // right           | PY |
        NY,     // bottom     +----+----+----+----+
        PY,     // top        | NX | PZ | PX | NZ |
        NZ,     // back       +----+----+----+----+
        PZ      // front           | NY |
                //                 +----+
    };

    struct Address
    {
        Face face;
        float s = 0;
        float t = 0;
    };

    Image faces[6];

    void setImageForFace(Face face, const Image& image);
    void getDirectionFor(float3& direction, Face face, int x, int y) const;
    void getDirectionFor(float3& direction, Face face, float x, float y) const;
    void getPixelFromDirection(float3& pixel, const float3& direction) const;
    void setAllFacesFromCross(const Image& image);
    Address getAddressFor(const float3& r) const;

    void makeSeamless();
};

inline void Cubemap::getDirectionFor(float3& direction, Face face, int x, int y) const
{
    getDirectionFor(direction, face, x + 0.5f, y + 0.5f);
}

inline void Cubemap::getDirectionFor(float3& direction, Face face, float x, float y) const
{
    const float scale = 2.0 / size;
    // map [0, dim] to [-1,1] with (-1,-1) at bottom left
    float cx = (x * scale) - 1;
    float cy = 1 - (y * scale);

    const float l = sqrt(cx * cx + cy * cy + 1);
    switch (face)
    {
    case PX:
        direction.set(1, cy, -cx);
        break;
    case NX:
        direction.set(-1, cy, cx);
        break;
    case PY:
        direction.set(cx, 1, -cy);
        break;
    case NY:
        direction.set(cx, -1, cy);
        break;
    case PZ:
        direction.set(cx, cy, 1);
        break;
    case NZ:
        direction.set(-cx, cy, -1);
        break;
    }

    direction *= 1 / l;
}
