#pragma once

#include "utils.h"
#include "Spherical.h"
#include "CubemapMipMap.h"
#include "Image.h"

namespace pkg {

enum ImageType { background = 0, brdf_ue4, thumbnail, specular_ue4, ImageTypeMax };
enum ImageEncoding { srgb = 0, luv, rgbm, rgbe, ImageEncodingMax };
enum ImageFormat { lut = 0, panorama, cubemap, ImageFormatMax };

struct ImageDescription
{
    Path filename;
    char name[64];
    int sizeUncompressed = 0;
    int sizeCompressed = 0;
    int width = 0;
    int height = 0;
    ImageEncoding encoding;
    ImageType type;
    ImageFormat format;
};

struct Package
{
    int version = 2;
    int numImages = 0;
    ImageDescription images[32];
    Spherical* diffuseSPH;
    const char* distDir;
    int needCompression = true;
    int padding;

    Package(const char* distDir) : distDir(distDir) {}

    void setSph(Spherical* sph) { diffuseSPH = sph; }

    // il faut ecrire directement et garder les infos filename ...
    void addPrefilterCubemap(const CubemapMipMap& cm, ImageEncoding encoding);
    void addPrefilterEquirectangular(const Image& img, ImageEncoding encoding);
    void addThumbnail(const Image& img);
    void addBackground(const Cubemap& cm, ImageEncoding encoding);
    // void addBackgroundCubemap(Cubemap* cm, const char* encoding)
    // {
    //     ImageDescription& image = images[numImages++];
    //     image.type = "background";
    //     image.encoding = encoding;
    //     image.format = "cubemap";
    //     image.width = image.height = cm.size;
    // }

    // void setLUT(Image* img)
    // {
    //     ImageDescription& image = images[numImages++];
    //     image.type = "brdf_ue4";
    //     image.encoding = "rg16";
    //     image.format = "lut";
    //     image.width = img.width;
    //     image.height = img.height;
    // }

    void write();
};

} // namespace pkg
