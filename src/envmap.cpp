#include "Cubemap.h"
#include "Image.h"
#include "Spherical.h"
#include "envUtils.h"
#include <getopt.h>
#include <stdio.h>
#include <thread>

bool debug = false;

int main(int argc, char** argv)
{
    char* file = argv[1];

    if (argc < 2)
    {
        printf("usage %s environment\n", argv[0]);
        return 1;
    }

    Image image;
    const char* distDir = "test";

    int loadImageResult = envUtils::loadImage(image, file);
    if (loadImageResult == 0)
    {

        int nbThread;
        nbThread = std::thread::hardware_concurrency();
        printf("using %d threads\n", nbThread);


        // needs to be sure that the panorama is pow of 2
        Cubemap cm;
        envUtils::createCubemap(cm, 256);

        envUtils::clampImage(image, 255);

        envUtils::writeThumbnail(distDir, "thumbnail", image, 256, 128);
        envUtils::equirectangularToCubemap(cm, image, nbThread);

        if (debug)
            envUtils::writeCubemap_hdr(distDir, "input", cm);

        CubemapMipMap cmMipMap;
        envUtils::createCubemapMipMap(cmMipMap, cm);

        if (debug)
            envUtils::writeCubemapMipMapFaces_hdr(distDir, "mipmap", cmMipMap);

        CubemapMipMap cmPrefilter;
        envUtils::prefilterCubemapGGX(cmPrefilter, cmMipMap, 4096, nbThread);

        if (debug)
            envUtils::writeCubemapMipMapFaces_hdr(distDir, "prefilter", cmPrefilter);

        Spherical spherical;
        envUtils::computeSphericalHarmonicsFromCubemap(spherical, cm);
        envUtils::writeSpherical_json(distDir, "spherical", spherical);
        envUtils::writeCubemapMipMap_luv(distDir, "prefilter", cmPrefilter);

        // CubemapMipMap cmDecode;
        // envUtils::readCubemapMipMap_luv(cmDecode, "test/prefilter.luv");
        // envUtils::writeCubemapMipMapFaces_hdr(distDir, "prefilter-decode", cmDecode);

        envUtils::writeCubemap_luv(distDir, "background", cmPrefilter.levels[2]);

        envUtils::freeImage(image);
        envUtils::freeCubemap(cm);
    }

    return loadImageResult;
}
