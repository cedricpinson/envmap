#pragma once

#include "float3.h"
#include "Image.h"
#include <math.h>

// major axis
// direction     target                              sc     tc    ma
// ----------    ---------------------------------   ---    ---   ---
//  +rx          GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT   -rz    -ry   rx
//  -rx          GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT   +rz    -ry   rx
//  +ry          GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT   +rx    +rz   ry
//  -ry          GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT   +rx    -rz   ry
//  +rz          GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT   +rx    -ry   rz
//  -rz          GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT   -rx    -ry   rz
// s   =   ( sc/|ma| + 1 ) / 2
// t   =   ( tc/|ma| + 1 ) / 2

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

    Image faces[6];

    void setImageForFace(Face face, const Image& image);
    double3 getDirectionFor(Face face, int x, int y) const;
    double3 getDirectionFor(Face face, double x, double y) const;
    void setAllFacesFromCross(const Image& image);

};

inline double3 Cubemap::getDirectionFor(Face face, int x, int y) const
{
    return getDirectionFor(face, x + 0.5, y + 0.5);
}

inline double3 Cubemap::getDirectionFor(Face face, double x, double y) const
{
    const double scale = 2.0 / size;
    // map [0, dim] to [-1,1] with (-1,-1) at bottom left
    double cx = (x * scale) - 1;
    double cy = 1 - (y * scale);

    double3 dir;
    const double l = sqrt(cx * cx + cy * cy + 1);
    switch (face)
    {
    case PX:
        dir = double3(1, cy, -cx);
        break;
    case NX:
        dir = double3(-1, cy, cx);
        break;
    case PY:
        dir = double3(cx, 1, -cy);
        break;
    case NY:
        dir = double3(cx, -1, cy);
        break;
    case PZ:
        dir = double3(cx, cy, 1);
        break;
    case NZ:
        dir = double3(-cx, cy, -1);
        break;
    }
    return dir * (1 / l);
}
