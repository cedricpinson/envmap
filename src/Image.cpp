#include "Image.h"

void Image::subset(const Image& image, int x, int y, int w, int h)
{
    rowInFloat3 = image.rowInFloat3;
    width = w;
    height = h;
    type = image.type;
    data = (float3*)&image.getPixel(x, y);
}



