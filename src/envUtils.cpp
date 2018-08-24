#include "envUtils.h"
#include "Cubemap.h"
#include "Spherical.h"
#include "float3.h"
#include "log.h"
#include "threadLines.h"
#include "utils.h"

#include <math.h>
#include <string.h>

// create directory
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr/tinyexr.h"

namespace envUtils {

int writeSpherical_json(const char* outputFilename, double* spherical)
{
    printf("writing file %s\n", outputFilename);

    const size_t BufferSize = 4096;
    char outputBuffer[BufferSize + 1];
    memset(outputBuffer, 0, BufferSize + 1);

    size_t count = 0;

    count += snprintf(&outputBuffer[count], BufferSize - count, "{\n");

#if 0
    // shRGB[I] * BandFactor[I]
    count += snprintf(&outputBuffer[count], BufferSize - count, "    \"shR\": [ %f", SHr[0] * SHBandFactor[0]);
    for (int i = 1; i < NUM_SH_COEFFICIENT; ++i)
        count += snprintf(&outputBuffer[count], BufferSize - count, ", %f", SHr[i] * SHBandFactor[i]);
    count += snprintf(&outputBuffer[count], BufferSize - count, "],\n");

    count += snprintf(&outputBuffer[count], BufferSize - count, "    \"shG\": [ %f", SHg[0] * SHBandFactor[0]);
    for (int i = 1; i < NUM_SH_COEFFICIENT; ++i)
        count += snprintf(&outputBuffer[count], BufferSize - count, ", %f", SHg[i] * SHBandFactor[i]);
    count += snprintf(&outputBuffer[count], BufferSize - count, " ],\n");

    count += snprintf(&outputBuffer[count], BufferSize - count, "    \"shB\": [ %f", SHb[0] * SHBandFactor[0]);
    for (int i = 0; i < NUM_SH_COEFFICIENT; ++i)
        count += snprintf(&outputBuffer[count], BufferSize - count, ", %f", SHb[i] * SHBandFactor[i]);
    count += snprintf(&outputBuffer[count], BufferSize - count, " ],\n\n");
#endif

    count += snprintf(&outputBuffer[count], BufferSize - count, "    \"shCoef\": [ %f, %f, %f", spherical[0],
                      spherical[1], spherical[2]);
    for (int i = 1; i < NUM_SH_COEFFICIENT; ++i)
    {
        count += snprintf(&outputBuffer[count], BufferSize - count, ", %f, %f,%f", spherical[i * 3],
                          spherical[i * 3 + 1], spherical[i * 3 + 2]);
    }
    count += snprintf(&outputBuffer[count], BufferSize - count, " ]\n");
    count += snprintf(&outputBuffer[count], BufferSize - count, "}\n");

    printf("%s\n", outputBuffer);
    FILE* fp = fopen(outputFilename, "w");
    if (!fp)
    {
        printf("can't write file %s", outputFilename);
        return 1;
    }
    if (fprintf(fp, "%s", outputBuffer) < 0)
    {
        printf("error writing file %s", outputFilename);
        return 1;
    }
    fclose(fp);

    return 0;
}

int load_image_exr(const char* filename, Image& image)
{
    const char* err;
    float* rgba;
    int ret = LoadEXR(&rgba, &image.width, &image.height, filename, &err);
    if (ret != 0)
    {
        printf("%s", err);
        return 1;
    }
    float3* data = (float3*)malloc(image.width * image.height * sizeof(float3));
    for (int i = 0; i < image.width * image.height; i++)
    {
        float3& texel = (float3&)data[i];
        texel[0] = rgba[i * 4];
        texel[1] = rgba[i * 4 + 1];
        texel[2] = rgba[i * 4 + 2];
    }
    image.data = data;
    image.rowInFloat3 = image.width;
    free(rgba);
    return 0;
}

void free_image_exr(Image& image) { free(image.data); }

int load_image_hdr(const char* filename, Image& image)
{

    int x, y, n;
    float* data = stbi_loadf(filename, &x, &y, &n, 3);
    if (!data)
        return 1;

    image.data = (float3*)data;
    image.width = x;
    image.height = y;
    image.rowInFloat3 = image.width;
    return 0;
}

void free_image_hdr(Image& image) { stbi_image_free(image.data); }

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

int writeImage_hdr(const char* filename, const Image& image)
{
    Path path;
    strncpy(path, filename, 1024);
    const char* extension = get_filename_ext(path);
    if (strlen(extension) == 0)
    {
        strcat(path, ".hdr");
    }

    printf("writing %s file\n", path);
    float3* data = image.data;
    if (image.width == image.rowInFloat3)
    {
        return stbi_write_hdr(path, image.width, image.height, 3, (float*)image.data);
    }

    // if the image is a subset of another image extract the subset to write it
    data = new float3[image.width * image.height];
    for (int y = 0; y < image.height; y++)
    {
        for (int x = 0; x < image.width; x++)
        {
            data[x + y * image.width] = image.data[x + y * image.rowInFloat3];
        }
    }
    int ret = stbi_write_hdr(path, image.width, image.height, 3, (float*)data);
    delete[] data;
    return ret;
}

int loadImage(Image& image, const char* filename)
{
    const char* ext = get_filename_ext(filename);
    printf("reading file %s\n", filename);
    if (strncmp(ext, "exr", 3) == 0)
    {
        image.type = Image::EXR;
        return load_image_exr(filename, image);
    }
    else if (strncmp(ext, "hdr", 3) == 0)
    {
        image.type = Image::HDR;
        return load_image_hdr(filename, image);
    }
    printf("file %s : extension not recognized, only hdr and exr files are supported", filename);
    return 0;
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

const char* getCubemapFaceName(Cubemap::Face face)
{
    switch (face)
    {
    case Cubemap::NX:
        return "face-NX";
    case Cubemap::PX:
        return "face-PX";
    case Cubemap::NY:
        return "face-NY";
    case Cubemap::PY:
        return "face-PY";
    case Cubemap::NZ:
        return "face-NZ";
    case Cubemap::PZ:
        return "face-PZ";
    }
    return "NONE";
}

void writeCubemap_hdr(const char* dir, const char* basename, const Cubemap& cm)
{
    make_directory(dir);

    Path path;
    char filename[512];
    for (int f = 0; f < 6; f++)
    {
        snprintf(filename, 511, "%s_%s", basename, envUtils::getCubemapFaceName((Cubemap::Face)f));
        create_path(path, dir, filename);
        const Image& image = cm.faces[f];
        writeImage_hdr(path, image);
    }

    snprintf(filename, 511, "%s_%s", basename, "cross");
    create_path(path, dir, filename);
    writeImage_hdr(path, cm.image);
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

struct EquirectangularProcessContext
{
    Cubemap* dst;
    Cubemap::Face face;
    const Image* src;
    EquirectangularProcessContext(Cubemap* dst, Cubemap::Face face, const Image* src)
        : dst(dst)
        , face(face)
        , src(src)
    {}
};

// this part is from filament https://github.com/google/filament/blob/master/tools/cmgen/src/CubemapUtils.cpp
void equirectangularToCubemapLines(EquirectangularProcessContext context, int startY, int stopY)
{
    Cubemap& dst = *context.dst;
    Cubemap::Face faceIdex = context.face;
    Image& face = dst.faces[faceIdex];
    const Image& src = *context.src;

    const size_t width = src.width;
    const size_t height = src.height;
    const double r = width * 0.5 * M_1_PI;
    int dim = face.width;

    for (int y = startY; y <= stopY; y++)
    {
        float3* data = &face.getPixel(0, y);
        for (int x = 0; x < dim; ++x, ++data)
        {
            // calculate how many samples we need based on dx, dy in the source
            // x =  cos(phi) sin(theta)
            // y = -sin(phi)
            // z = -cos(phi) cos(theta)
            float3 s0;
            dst.getDirectionFor(s0, faceIdex, x, y);
            const double t0 = atan2(s0[0], -s0[2]);
            const double p0 = asin(s0[1]);
            float3 s1;
            dst.getDirectionFor(s1, faceIdex, x + 1, y + 1);
            const double t1 = atan2(s1[0], -s1[2]);
            const double p1 = asin(s1[1]);
            const double dt = abs(t1 - t0);
            const double dp = abs(p1 - p0);
            const double dx = abs(r * dt);
            const double dy = abs(r * dp * s0[1]);
            const size_t numSamples = (size_t const)ceil(fmax(dx, dy));
            const float iNumSamples = 1.0f / numSamples;

            double3 c = double3(0, 0, 0);
            for (size_t sample = 0; sample < numSamples; sample++)
            {
                // Generate numSamples in our destination pixels and map them to input pixels
                const double2 h = hammersley(size_t(sample), iNumSamples);
                float3 s;
                dst.getDirectionFor(s, faceIdex, float(x + h[0]), float(y + h[1]));
                double xf = atan2(s[0], -s[2]) * M_1_PI; // range [-1.0, 1.0]
                double yf = asin(-s[1]) * (2 * M_1_PI);  // range [-1.0, 1.0]
                xf = (xf + 1) * 0.5 * (width - 1);       // range [0, width [
                yf = (yf + 1) * 0.5 * (height - 1);      // range [0, height[
                // we can't use filterAt() here because it reads past the width/height
                // which is okay for cubmaps but not for square images
                const float3 pixel = src.getPixel((uint32_t)xf, (uint32_t)yf);
                c[0] += pixel[0];
                c[1] += pixel[1];
                c[2] += pixel[2];
            }
            c *= iNumSamples;
            float3& resultPixel = *data;
            resultPixel[0] = c[0];
            resultPixel[1] = c[1];
            resultPixel[2] = c[2];
        }
    }
}

void equirectangularToCubemap(Cubemap& dst, const Image& src)
{
    int nbThread;
    nbThread = std::thread::hardware_concurrency();
    printf("using %d threads\n", nbThread);
    auto t = logStart("equirectangularToCubemap");

    for (int f = 0; f < 6; f++)
    {
        EquirectangularProcessContext context(&dst, (Cubemap::Face)f, &src);
        threadLines(equirectangularToCubemapLines, context, dst.size, nbThread);
    }
    logEnd(t);
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
                srcFace.filterAt(dstLine[x], x * scale + 0.5, y * scale + 0.5);
            }
        }
    }
}

void createCubemapMipMap(CubemapMipMap& cmMipMap, const Cubemap& cm)
{
    size_t maxSize = cm.size;
    size_t numMipMap = log2(maxSize) + 1;

    cmMipMap.init(numMipMap);

    // copy first level 0
    envUtils::createCubemap(cmMipMap.levels[0], maxSize);
    memcpy(cmMipMap.levels[0].image.data, cm.image.data, cm.image.width * cm.image.height * sizeof(float3));
    cmMipMap.levels[0].makeSeamless();

    for (size_t i = 1; i < numMipMap; i++)
    {
        size_t size = pow(2, (numMipMap - 1) - i);
        envUtils::createCubemap(cmMipMap.levels[i], size);
        envUtils::downsampleCubemapLevelBoxFilter(cmMipMap.levels[i], cmMipMap.levels[i - 1]);
        cmMipMap.levels[i].makeSeamless();
    }
}

int writeCubemapMipMap_hdr(const char* dir, const char* basename, const CubemapMipMap& cm)
{
    make_directory(dir);

    Path path;
    char filename[256];
    for (int i = 0; i < cm.numLevel; i++)
    {
        snprintf(filename, 255, "%s_level_%d", basename, i);
        create_path(path, dir, filename);
        const Image& image = cm.levels[i].image;
        writeImage_hdr(path, image);
    }
    return 0;
}

int writeCubemapMipMapFaces_hdr(const char* dir, const char* basename, const CubemapMipMap& cmMipMap)
{
    make_directory(dir);

    Path path;
    char filename[256];
    for (int i = 0; i < cmMipMap.numLevel; i++)
    {
        snprintf(filename, 255, "%s_level_%d", basename, i);
        create_path(path, dir, filename);
        const Cubemap& cm = cmMipMap.levels[i];
        writeCubemap_hdr(path, filename, cm);
    }
    return 0;
}

} // namespace envUtils
