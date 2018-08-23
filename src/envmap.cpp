#include "Cubemap.h"
#include "Image.h"
#include "Spherical.h"
#include "envUtils.h"
#include <getopt.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    char* file = argv[1];

    if (argc < 2)
    {
        printf("usage %s environment\n", argv[0]);
        return 1;
    }

    Image image;

    if (envUtils::loadImage(image, file) == 0)
    {

        Cubemap cm;
        envUtils::createCubemap(cm, 256);

        envUtils::clampImage(image, 255);
        envUtils::equirectangularToCubemap(cm, image);
        envUtils::writeCubemap_hdr("test", "input", cm);

        CubemapMipMap cmMipMap;
        envUtils::createCubemapMipMap(cmMipMap, cm);
        envUtils::writeCubemapMipMapFaces_hdr("test", "mipmap", cmMipMap);

        CubemapMipMap cmPrefilter;
        envUtils::prefilterCubemapGGX(cmPrefilter, cmMipMap, 1024);
        envUtils::writeCubemapMipMapFaces_hdr("test", "prefilter", cmPrefilter);

        envUtils::writeImage_hdr("./test.hdr", image);

        Spherical spherical;
        envUtils::computeSphericalHarmonicsFromCubemap(spherical, cm);
        envUtils::writeSpherical_json("./test/spherical.json", spherical);

        envUtils::freeImage(image);
        envUtils::freeCubemap(cm);
    }

    return 0;
}
