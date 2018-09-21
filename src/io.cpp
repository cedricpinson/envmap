#include "io.h"
#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "Spherical.h"
#include "encode.h"
#include "envUtils.h"
#include "utils.h"

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

namespace io {

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

int writeSpherical_json(const char* path, double* spherical)
{
    printf("writing file %s\n", path);

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
    FILE* fp = fopen(path, "w");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }
    if (fprintf(fp, "%s", outputBuffer) < 0)
    {
        printf("error writing file %s", path);
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

int writeImage_hdr(const char* path, const Image& image)
{
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

int writeImage_luv(const char* path, const Image& image)
{
    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    printf("writing %s file\n", path);

    int w = image.width;
    int h = image.height;
    uint8_t* dest = new uint8_t[w * h * 4];
    uint8_t rgba8[4];
    for (int y = 0; y < image.height; y++)
    {
        const float3* lineSrc = &image.getPixel(0, y);
        for (int x = 0; x < image.width; x++)
        {
            // interleave the face
            encodeLUV(rgba8, lineSrc[x].ptr());
            dest[y * w + x + (w * h * 0)] = rgba8[0];
            dest[y * w + x + (w * h * 1)] = rgba8[1];
            dest[y * w + x + (w * h * 2)] = rgba8[2];
            dest[y * w + x + (w * h * 3)] = rgba8[3];
        }
    }

    fwrite(dest, w * h * 4, 1, fp);
    delete[] dest;
    return 0;
}

int writeImage_rgbm(const char* path, const Image& image)
{
    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    printf("writing %s file\n", path);

    int w = image.width;
    int h = image.height;
    uint8_t* dest = new uint8_t[w * h * 4];
    for (int y = 0; y < image.height; y++)
    {
        const float3* lineSrc = &image.getPixel(0, y);
        uint8_t* lineDst = &dest[y * w * 4];
        for (int x = 0; x < image.width; x++, lineDst += 4)
        {
            encodeRGBM(lineDst, lineSrc[x].ptr());
        }
    }

    fwrite(dest, w * h * 4, 1, fp);
    delete[] dest;
    return 0;
}

int writeImage_ldr(const char* path, const Image& image)
{
    printf("writing %s file\n", path);
    uint8_t* dest = new uint8_t[image.width * image.height * 3];

    for (int y = 0; y < image.height; y++)
    {
        for (int x = 0; x < image.width; x++)
        {
            uint8_t* pixelDst = &dest[(x + y * image.width) * 3];
            const float3& pixelSrc = image.getPixel(x, y);
            // operation hdr to ldr
            tonemap(pixelDst, pixelSrc.ptr());
        }
    }
    int ret = stbi_write_jpg(path, image.width, image.height, 3, dest, 92);

    delete[] dest;
    return ret;
}

int loadImage(Image& image, const char* filename)
{
    const char* ext = getFilenameExtension(filename);
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

#define ENTERLEAVE
int readCubemapMipMap_luv(CubemapMipMap& cmMipMap, const char* path)
{
    int baseSize = 256;
    int numLod = (int)log2(baseSize) + 1;

    cmMipMap.init(numLod);

    FILE* fp = fopen(path, "rb");
    if (!fp)
    {
        printf("can't open file %s", path);
        return 1;
    }

    uint8_t* src = new uint8_t[baseSize * baseSize * 4];
    for (int mipLevel = 0; mipLevel < numLod; mipLevel++)
    {
        int size = (int)pow(2, (numLod - mipLevel - 1));
        Cubemap cm;
        envUtils::createCubemap(cm, size);

        cmMipMap.levels[mipLevel] = cm;
        for (int f = 0; f < 6; f++)
        {
            fread(src, size * size * 4, 1, fp);
            for (int y = 0; y < size; y++)
            {
                float3* lineDst = &cm.faces[f].getPixel(0, y);
                for (int x = 0; x < size; x++)
                {
#ifdef ENTERLEAVE
                    uint8_t rgba8[4];
                    // interleave the face
                    rgba8[0] = src[y * size + x + (size * size * 0)];
                    rgba8[1] = src[y * size + x + (size * size * 1)];
                    rgba8[2] = src[y * size + x + (size * size * 2)];
                    rgba8[3] = src[y * size + x + (size * size * 3)];
                    decodeLUV(lineDst[x].ptr(), rgba8);
#else
                    decodeLUV(lineDst[x].ptr(), &src[(y * size + x) * 4]);
#endif
                }
            }
        }
    }
    fclose(fp);
    delete[] src;
    return 0;
}

static void writeCubemap_luv_internal(FILE* fp, uint8_t* tempBuffer, const Cubemap& cm)
{
    uint8_t* dest = tempBuffer;
    int size = cm.size;
    for (int i = 0; i < 6; i++)
    {
        for (int y = 0; y < size; y++)
        {
            const float3* lineSrc = &cm.faces[i].getPixel(0, y);
            for (int x = 0; x < size; x++)
            {
#ifdef ENTERLEAVE
                uint8_t rgba8[4];
                // interleave the face
                encodeLUV(rgba8, lineSrc[x].ptr());
                dest[y * size + x + (size * size * 0)] = rgba8[0];
                dest[y * size + x + (size * size * 1)] = rgba8[1];
                dest[y * size + x + (size * size * 2)] = rgba8[2];
                dest[y * size + x + (size * size * 3)] = rgba8[3];
#else
                encodeLUV(&Dst[(y * size + x) * 4], lineSrc[x].ptr());
#endif
            }
        }
        fwrite(dest, size * size * 4, 1, fp);
    }
}

int writeCubemap_luv(const char* path, const Cubemap& cm)
{
    printf("writing %s file\n", path);

    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    // allocate only one face
    int maxSize = cm.size;
    uint8_t* dest = new uint8_t[maxSize * maxSize * 4];

    writeCubemap_luv_internal(fp, dest, cm);

    fclose(fp);
    delete[] dest;
    return 0;
}

static void writeCubemap_rgbm_internal(FILE* fp, uint8_t* tempBuffer, const Cubemap& cm)
{
    uint8_t* dest = tempBuffer;
    int size = cm.size;
    for (int i = 0; i < 6; i++)
    {
        for (int y = 0; y < size; y++)
        {
            const float3* lineSrc = &cm.faces[i].getPixel(0, y);
            for (int x = 0; x < size; x++)
            {
                encodeRGBM(&dest[(y * size + x) * 4], lineSrc[x].ptr());
            }
        }
        fwrite(dest, size * size * 4, 1, fp);
    }
}

int writeCubemap_rgbm(const char* path, const Cubemap& cm)
{
    printf("writing %s file\n", path);

    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    // allocate only one face
    int maxSize = cm.size;
    uint8_t* dest = new uint8_t[maxSize * maxSize * 4];

    writeCubemap_rgbm_internal(fp, dest, cm);

    fclose(fp);
    delete[] dest;
    return 0;
}

int writeCubemapMipMap_rgbm(const char* path, const CubemapMipMap& cmMipMap)
{
    printf("writing %s file\n", path);

    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    // allocate only one face
    int maxSize = cmMipMap.levels[0].size;
    uint8_t* dest = new uint8_t[maxSize * maxSize * 4];

    for (int mipLevel = 0; mipLevel < cmMipMap.numLevel; mipLevel++)
    {
        const Cubemap& cm = cmMipMap.levels[mipLevel];
        writeCubemap_rgbm_internal(fp, dest, cm);
    }
    fclose(fp);
    delete[] dest;
    return 0;
}

int writeCubemapMipMap_luv(const char* path, const CubemapMipMap& cmMipMap)
{
    printf("writing %s file\n", path);

    FILE* fp = fopen(path, "wb");
    if (!fp)
    {
        printf("can't write file %s", path);
        return 1;
    }

    // allocate only one face
    int maxSize = cmMipMap.levels[0].size;
    uint8_t* dest = new uint8_t[maxSize * maxSize * 4];

    for (int mipLevel = 0; mipLevel < cmMipMap.numLevel; mipLevel++)
    {
        const Cubemap& cm = cmMipMap.levels[mipLevel];
        writeCubemap_luv_internal(fp, dest, cm);
    }
    fclose(fp);
    delete[] dest;
    return 0;
}

int writeCubemapMipMapDebug(const char* dir, const char* basename, const CubemapMipMap& cmMipMap)
{
    makeDirectory(dir);

    Path path;
    char filename[256];
    for (int i = 0; i < cmMipMap.numLevel; i++)
    {
        snprintf(filename, 255, "%s_level_%d", basename, i);
        createPath(path, dir, filename);
        const Cubemap& cm = cmMipMap.levels[i];
        writeCubemapDebug(path, filename, cm);
    }
    return 0;
}

void writeCubemapDebug(const char* dir, const char* basename, const Cubemap& cm)
{
    makeDirectory(dir);

    char filename[512];
    for (int f = 0; f < 6; f++)
    {
        snprintf(filename, 511, "%s/%s_%s", dir, basename, getCubemapFaceName((Cubemap::Face)f));
        const Image& image = cm.faces[f];
        writeImage_hdr(filename, image);
    }

    snprintf(filename, 511, "%s/%s_%s", dir, basename, "cross");
    writeImage_hdr(filename, cm.image);
}

} // namespace io
