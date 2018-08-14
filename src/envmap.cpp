#include "Cubemap.h"
#include "Image.h"
#include "Spherical.h"
#include "envUtils.h"
#include <getopt.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    char* file = argv[1];

    Image image;

    if (envUtils::loadImage(image, file) == 0)
    {

        Cubemap cm;
        envUtils::createCubemap(cm, 256);

        envUtils::equirectangularToCubemap(cm, image);
        envUtils::writeCubemap_hdr("test", cm);

        CubemapMipMap cmMipMap;
        envUtils::createCubemapMipMap(cmMipMap, cm);
        envUtils::writeCubemapMipMap_hdr("test", cmMipMap);

        CubemapMipMap cmPrefilter;
        envUtils::prefilterCubemapGGX(cmPrefilter, cmMipMap, 512);

        envUtils::writeImage_hdr("./test.hdr", image);

        Spherical spherical;
        envUtils::computeSphericalHarmonicsFromCubemap(spherical, cm);
        envUtils::writeSpherical_json("./test/spherical.json", spherical);

        envUtils::freeImage(image);
        envUtils::freeCubemap(cm);
    }

    return 0;
}
