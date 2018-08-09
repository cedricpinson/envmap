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
