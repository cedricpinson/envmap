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

    void getDirectionFixUpFor(float* direction, Face face, int x, int y) const;

    void setAllFacesFromCross(const Image& image);

    void makeSeamless();

    void _remapUV(float& cx, float& cy, float x, float y) const;
    void _remapFixUpUV(float& cx, float& cy, float x, float y) const;
    void _getDirectionFor(float* direction, Face face, float cx, float cy) const;

    inline static void getAddressFor(Address& addr, const float* r)
    {
        float sc, tc, ma;
        const float rx = fabsf(r[0]);
        const float ry = fabsf(r[1]);
        const float rz = fabsf(r[2]);
        if (rx >= ry && rx >= rz)
        {
            ma = rx;
            if (r[0] >= 0)
            {
                addr.face = Cubemap::Face::PX;
                sc = -r[2];
                tc = -r[1];
            }
            else
            {
                addr.face = Cubemap::Face::NX;
                sc = r[2];
                tc = -r[1];
            }
        }
        else if (ry >= rx && ry >= rz)
        {
            ma = ry;
            if (r[1] >= 0)
            {
                addr.face = Cubemap::Face::PY;
                sc = r[0];
                tc = r[2];
            }
            else
            {
                addr.face = Cubemap::Face::NY;
                sc = r[0];
                tc = -r[2];
            }
        }
        else
        {
            ma = rz;
            if (r[2] >= 0)
            {
                addr.face = Cubemap::Face::PZ;
                sc = r[0];
                tc = -r[1];
            }
            else
            {
                addr.face = Cubemap::Face::NZ;
                sc = -r[0];
                tc = -r[1];
            }
        }
        // ma is guaranteed to be >= sc and tc
        addr.s = (sc / ma + 1.f) * 0.5f;
        addr.t = (tc / ma + 1.f) * 0.5f;
    }
};

inline void Cubemap::getDirectionFor(float* direction, Face face, int x, int y) const
{
    getDirectionFor(direction, face, x + 0.5f, y + 0.5f);
}

inline void Cubemap::_remapUV(float& cx, float& cy, float x, float y) const
{
    const float scale = 2.0f / size;
    // map [0, dim] to [-1,1] with (-1,-1) at bottom left
    cx = (x * scale) - 1;
    cy = 1 - (y * scale);
}

inline void Cubemap::_remapFixUpUV(float& cx, float& cy, float x, float y) const
{
    if (size == 1)
    {
        cx = 0.0f;
        cy = 0.0f;
        return;
    }
    // transform from [0..res - 1] to [-1 .. 1], match up edges exactly.
    float scale = 2.0f / (size - 1.0f);
    cx = (x * scale) - 1;
    cy = 1 - (y * scale);
}

inline void Cubemap::_getDirectionFor(float* direction, Face face, float cx, float cy) const
{
    const float l = 1.0f / sqrtf(cx * cx + cy * cy + 1.0f);
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

    direction[0] *= l;
    direction[1] *= l;
    direction[2] *= l;
}

inline void Cubemap::getDirectionFor(float* direction, Face face, float x, float y) const
{
    float cx, cy;
    _remapUV(cx, cy, x, y);
    _getDirectionFor(direction, face, cx, cy);
}

inline void Cubemap::getDirectionFixUpFor(float* direction, Face face, int x, int y) const
{
    float cx, cy;
    _remapFixUpUV(cx, cy, (float)x, (float)y);
    _getDirectionFor(direction, face, cx, cy);
}
