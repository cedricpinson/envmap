#include "Cubemap.h"

void setFaceFromCross(Cubemap& cm, Cubemap::Face face, const Image& image)
{
    int dim = cm.size;
    int x = 0;
    int y = 0;

    switch (face)
    {
    case Cubemap::NX:
        x = 0, y = dim;
        break;
    case Cubemap::PX:
        x = 2 * dim, y = dim;
        break;
    case Cubemap::NY:
        x = dim, y = 2 * dim;
        break;
    case Cubemap::PY:
        x = dim, y = 0;
        break;
    case Cubemap::NZ:
        x = 3 * dim, y = dim;
        break;
    case Cubemap::PZ:
        x = dim, y = dim;
        break;
    }

    Image subImage;
    subImage.subset(image, x, y, dim, dim);
    cm.setImageForFace(face, subImage);
}

void Cubemap::setImageForFace(Face face, const Image& image) { faces[face] = image; }

void Cubemap::setAllFacesFromCross(const Image& image)
{
    setFaceFromCross(*this, Cubemap::NX, image);
    setFaceFromCross(*this, Cubemap::PX, image);
    setFaceFromCross(*this, Cubemap::NY, image);
    setFaceFromCross(*this, Cubemap::PY, image);
    setFaceFromCross(*this, Cubemap::NZ, image);
    setFaceFromCross(*this, Cubemap::PZ, image);
}

// From filament https://github.com/google/filament/blob/master/tools/cmgen/src/Cubemap.cpp#L97
/*
 * We handle "seamless" cubemaps by duplicating a row to the bottom, or column to the right
 * of each faces that don't have an adjacent face in the image (the duplicate is taken from the
 * adjacent face in the cubemap).
 * This is because when accessing an image with bilinear filtering, we always overshoot to the
 * right or bottom. This works well with cubemaps stored as a cross in memory.
 */
void Cubemap::makeSeamless()
{
    size_t dim = size;

    auto stich = [&](float3* dst, int incDst, float3 const* src, int incSrc) {
        for (size_t i = 0; i < dim; ++i)
        {
            *dst = *src;
            dst = dst + incDst;
            src = src + incSrc;
        }
    };

    const size_t bpr = faces[Face::NX].rowInFloat3;

    stich(&faces[Face::NX].getPixel(0, dim), 1, &faces[Face::NY].getPixel(0, dim - 1), -bpr);
    stich(&faces[Face::PY].getPixel(dim, 0), bpr, &faces[Face::PX].getPixel(dim - 1, 0), -1);

    stich(&faces[Face::PX].getPixel(0, dim), 1, &faces[Face::NY].getPixel(dim - 1, 0), bpr);

    stich(&faces[Face::NY].getPixel(dim, 0), bpr, &faces[Face::PX].getPixel(0, dim - 1), 1);

    // horizontal cross
    stich(&faces[Face::NZ].getPixel(0, dim), 1, &faces[Face::NY].getPixel(dim - 1, dim - 1), -1);

    stich(&faces[Face::NZ].getPixel(dim, 0), bpr, &faces[Face::NX].getPixel(0, 0), bpr);

    stich(&faces[Face::NY].getPixel(0, dim), 1, &faces[Face::NZ].getPixel(dim - 1, dim - 1), -1);
}
