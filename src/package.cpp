#include "package.h"
#include "envUtils.h"
#include <stdio.h>
#include <string.h>

namespace pkg {

static const char* getImageTypeStr(ImageType type)
{
    switch (type)
    {
    case background:
        return "background";
    case brdf_ue4:
        return "brdf_ue4";
    case thumbnail:
        return "thumbnail";
    case specular_ue4:
        return "specular_ue4";
    default:
        return "unknown";
    }
    return "unknown";
}

static const char* getImageEncodingStr(ImageEncoding encoding)
{
    switch (encoding)
    {
    case srgb:
        return "srgb";
    case luv:
        return "luv";
    case rgbm:
        return "rgbm";
    case rgbe:
        return "rgbe";
    default:
        return "unknown";
    }
    return "unknown";
}

static const char* getImageFormatStr(ImageFormat format)
{
    switch (format)
    {
    case lut:
        return "lut";
    case panorama:
        return "panorama";
    case cubemap:
        return "cubemap";
    default:
        return "unknown";
    }
    return "unknown";
}

int printImageHeader(char* buffer, const ImageDescription& image)
{
    int size = 0;
    size += sprintf(buffer + size, "      \"format\": \"%s\",\n", getImageFormatStr(image.format));
    size += sprintf(buffer + size, "      \"type\": \"%s\",\n", getImageTypeStr(image.type));
    size += sprintf(buffer + size, "      \"encoding\": \"%s\",\n", getImageEncodingStr(image.encoding));
    return size;
}

int printImage(char* buffer, const ImageDescription& image)
{
    int size = 0;
    size += sprintf(buffer + size, "        {\n");
    size += sprintf(buffer + size, "          \"name\": \"%s\",\n", image.name);
    size += sprintf(buffer + size, "          \"filename\": \"%s\",\n", image.filename);
    if (image.sizeCompressed)
    {
        size += sprintf(buffer + size, "          \"sizeCompressed\": %d,\n", image.sizeCompressed);
        size += sprintf(buffer + size, "          \"sizeUncompressed\": %d,\n", image.sizeUncompressed);
    }
    size += sprintf(buffer + size, "          \"width\": %d,\n", image.width);
    size += sprintf(buffer + size, "          \"height\": %d\n", image.height);
    size += sprintf(buffer + size, "        }\n");
    return size;
}

void Package::addPrefilterCubemap(const CubemapMipMap& cm, ImageEncoding encoding)
{
    Path filename;
    sprintf(filename, "specular_cubemap_%d_%s.bin", cm.levels[0].size, getImageEncodingStr(encoding));

    ImageDescription& imageDescription = images[numImages++];

    if (encoding == ImageEncoding::luv)
    {
        envUtils::writeCubemapMipMap_luv(distDir, filename, cm);
    }
    Path path;
    createPath(path, distDir, filename);

    int sizeUncompressed = getFileSize(path);
    int sizeCompressed = 0;

    memcpy(imageDescription.filename, filename, sizeof(filename));

    if (needCompression)
    {
        sizeCompressed = compressGZ(path);
        strcat(imageDescription.filename, ".gz");
    }

    imageDescription.type = ImageType::specular_ue4;
    imageDescription.encoding = encoding;
    imageDescription.format = ImageFormat::cubemap;
    imageDescription.width = imageDescription.height = cm.levels[0].size;
    imageDescription.sizeUncompressed = sizeUncompressed;
    imageDescription.sizeCompressed = sizeCompressed;
}

void Package::addPrefilterEquirectangular(const Image& image, ImageEncoding encoding)
{
    Path filename;
    sprintf(filename, "specular_panorama_%d_%s.bin", image.width, getImageEncodingStr(encoding));

    ImageDescription& imageDescription = images[numImages++];

    if (encoding == ImageEncoding::luv)
    {
        envUtils::writeImage_luv(distDir, filename, image);
    }

    Path path;
    createPath(path, distDir, filename);

    int sizeUncompressed = getFileSize(path);
    int sizeCompressed = 0;

    memcpy(imageDescription.filename, filename, sizeof(filename));

    if (needCompression)
    {
        sizeCompressed = compressGZ(path);
        strcat(imageDescription.filename, ".gz");
    }

    imageDescription.type = ImageType::specular_ue4;
    imageDescription.encoding = encoding;
    imageDescription.format = ImageFormat::panorama;
    imageDescription.width = image.width;
    imageDescription.height = image.height;
    imageDescription.sizeUncompressed = sizeUncompressed;
    imageDescription.sizeCompressed = sizeCompressed;
}

void Package::addThumbnail(const Image& image)
{
    Path filename;
    sprintf(filename, "thumbnail_%d.jpg", image.width);

    Path path;
    createPath(path, distDir, filename);

    ImageDescription& imageDescription = images[numImages++];

    envUtils::writeImage_ldr(distDir, filename, image);

    int sizeUncompressed = getFileSize(path);
    memcpy(imageDescription.filename, filename, sizeof(filename));

    imageDescription.type = ImageType::thumbnail;
    imageDescription.encoding = ImageEncoding::srgb;
    imageDescription.format = ImageFormat::panorama;
    imageDescription.width = image.width;
    imageDescription.height = image.height;
    imageDescription.sizeUncompressed = sizeUncompressed;
}

void Package::addBackground(const Cubemap& cm, ImageEncoding encoding)
{
    Path filename;
    sprintf(filename, "background_cubemap_%d_%s.bin", cm.size, getImageEncodingStr(encoding));

    ImageDescription& imageDescription = images[numImages++];

    if (encoding == ImageEncoding::luv)
    {
        envUtils::writeCubemap_luv(distDir, filename, cm);
    }

    Path path;
    createPath(path, distDir, filename);

    int sizeUncompressed = getFileSize(path);
    int sizeCompressed = 0;

    memcpy(imageDescription.filename, filename, sizeof(filename));

    if (needCompression)
    {
        sizeCompressed = compressGZ(path);
        strcat(imageDescription.filename, ".gz");
    }

    imageDescription.type = ImageType::background;
    imageDescription.encoding = encoding;
    imageDescription.format = ImageFormat::cubemap;
    imageDescription.width = imageDescription.height = cm.size;
    imageDescription.sizeUncompressed = sizeUncompressed;
    imageDescription.sizeCompressed = sizeCompressed;
}

void Package::write()
{

    // 64k config file
    char* config = new char[65536];

    int size = 0;
    size += sprintf(config, "{\n  \"version\": %d,\n", version);
    size += sprintf(config + size, "  \"textures\": [\n");

    bool firstImage = true;
    for (int type = 0; type < ImageTypeMax; type++)
    {
        for (int encoding = 0; encoding < ImageEncodingMax; encoding++)
        {
            for (int format = 0; format < ImageFormatMax; format++)
            {

                int nbImages = 0;

                for (int i = 0; i < numImages; i++)
                {
                    const ImageDescription& image = images[i];
                    if (image.format == (ImageFormat)format && image.encoding == (ImageEncoding)encoding &&
                        image.type == (ImageType)type)
                    {
                        if (!nbImages)
                        {
                            if (!firstImage)
                                size += sprintf(config + size, ",\n");
                            size += sprintf(config + size, "    {\n");
                            size += printImageHeader(config + size, image);
                            size += sprintf(config + size, "      \"images\": [\n");
                            firstImage = false;
                        }
                        size += printImage(config + size, image);
                        nbImages++;
                    }
                }
                if (nbImages)
                {
                    size += sprintf(config + size, "      ]\n");
                    size += sprintf(config + size, "    }");
                }
            }
        }
    }

    size += sprintf(config + size, "\n  ]\n");
    size += sprintf(config + size, "}\n");

    printf("%s", config);

    delete[] config;
}

} // namespace pkg
