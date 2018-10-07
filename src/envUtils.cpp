#include "envUtils.h"
#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "Spherical.h"
#include "float3.h"
#include "log.h"
#include "threadLines.h"

#include <math.h>
#include <string.h>

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#define STBIR_DEFAULT_FILTER_DOWNSAMPLE STBIR_FILTER_CATMULLROM
#include "stb/stb_image_resize.h"

namespace envUtils {

void free_image_exr(Image& image) { free(image.data); }
void free_image_hdr(Image& image) { free(image.data); }

void createImage(Image& image, int width, int height)
{
    size_t size = width * height * sizeof(float3);
    image.data = (float3*)malloc(size);
    memset(image.data, 0, size);
    image.width = width;
    image.rowInFloat3 = width;
    image.height = height;
    image.type = Image::RAW;
}

void resizeImage(Image& dst, const Image& image, int width, int height)
{
    auto t = logStart("resizeImage");

    envUtils::createImage(dst, width, height);
    stbir_resize_float(image.getData(), image.width, image.height, 0, dst.getData(), width, height, 0, 3);
    logEnd(t);
}

void freeImage(Image& image)
{
    if (image.type == Image::EXR)
    {
        free_image_exr(image);
    }
    else if (image.type == Image::HDR)
    {
        free_image_hdr(image);
    }
    else
    {
        ::free(image.data);
    }
}

// cubemap

void initCubemap(Cubemap& cm, const Image& image) { cm.image = image; }

void createCubemap(Cubemap& cm, int size)
{
    // use same technique as filament for bilinear filtering
    // https://github.com/google/filament/blob/master/tools/cmgen/src/Cubemap.cpp#L90
    int width = 4 * size + 1;
    int height = 3 * size + 1;
    createImage(cm.image, width, height);
    cm.size = size;
    cm.setAllFacesFromCross(cm.image);
}

void freeCubemap(Cubemap& cm)
{
    if (cm.image.data)
    {
        freeImage(cm.image);
    }
}

void clampImage(Image& src, float maxValue)
{
    const int size = src.width * src.height;
    float3* data = src.data;
    for (int i = 0; i < size; i++, data++)
    {
        float3& pixel = *data;
        pixel[0] = fmin(pixel[0], maxValue);
        pixel[1] = fmin(pixel[1], maxValue);
        pixel[2] = fmin(pixel[2], maxValue);
    }
}

void downsampleCubemapLevelBoxFilter(Cubemap& dst, const Cubemap& src)
{
    int scale = 2;
    int dim = dst.size;
    for (int f = 0; f < 6; f++)
    {
        const Image& srcFace = src.faces[f];

        for (int y = 0; y < dim; y++)
        {
            float3* dstLine = &dst.faces[f].getPixel(0, y);
            for (int x = 0; x < dim; ++x)
            {
                srcFace.filterAt(dstLine[x].ptr(), x * scale + 0.5, y * scale + 0.5);
            }
        }
    }
}

void createCubemapMipMap(CubemapMipMap& cmMipMap, const Cubemap& cm)
{
    auto t = logStart("createCubemapMipMap");
    size_t maxSize = cm.size;
    size_t numMipMap = (size_t)log2(maxSize) + 1;

    cmMipMap.init(numMipMap);

    // copy first level 0
    envUtils::createCubemap(cmMipMap.levels[0], maxSize);
    memcpy(cmMipMap.levels[0].image.data, cm.image.data, cm.image.width * cm.image.height * sizeof(float3));
    cmMipMap.levels[0].makeSeamless();

    for (size_t i = 1; i < numMipMap; i++)
    {
        size_t size = (size_t)pow(2, (numMipMap - 1) - i);
        envUtils::createCubemap(cmMipMap.levels[i], size);
        envUtils::downsampleCubemapLevelBoxFilter(cmMipMap.levels[i], cmMipMap.levels[i - 1]);
        cmMipMap.levels[i].makeSeamless();
    }
    logEnd(t);
}

void freeCubemapMipMap(CubemapMipMap& cmMipMap)
{
    for (int i = 0; i < cmMipMap.numLevel; i++)
    {
        freeCubemap(cmMipMap.levels[i]);
    }
    cmMipMap.numLevel = 0;
}

void packPrefilterCubemapToEquilateral(Image& equirectangular, const CubemapMipMap& src, int nbThreads)
{

    auto t = logStart("cubemapToEquirectangular");
    int numMipMap = src.numLevel;
    int cubemapSize = src.levels[0].size;

    // image that will contains mipmap of panorama
    envUtils::createImage(equirectangular, cubemapSize * 4, cubemapSize * 4);

    int sizeLevel;
    int yOffset = 0;
    Image equirectangularLevel;

    for (int i = 0; i < numMipMap; i++)
    {
        sizeLevel = src.levels[i].size;
        equirectangularLevel.subset(equirectangular, 0, yOffset, sizeLevel * 4, sizeLevel * 2);
        envUtils::cubemapToEquirectangular(equirectangularLevel, src.levels[i], nbThreads);
        yOffset += sizeLevel * 2;
    }

    logEnd(t);
}

} // namespace envUtils
