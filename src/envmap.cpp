#include "Cubemap.h"
#include "Image.h"
#include "Spherical.h"
#include "envUtils.h"
#include "utils.h"
#include <getopt.h>
#include <stdio.h>
#include <thread>

#include <getopt.h>

struct Options
{
    int numSamples = 1024;
    int cubemapSize = 256;
    const char* inputFilename;
    const char* distDir = "dist";
    int thumbnailSize = 256;
    int nbThreads = 0;
    bool debug = false;
    bool quiet = false;
    unsigned short padding;
    unsigned int padding2;
};

void printUsage()
{
    const char* text =
        "envmap is a tool to precompute environment for PBR rendering "
        "usages:"
        " envmap [options] panorama_env.hdr"
        ""
        "example"
        " envmap -s 512 -n 4096 panorama_env.hdr"
        ""
        "Options:"
        " -s [--ibl-size] to define the cubemap resolution computed, default 256"
        " -n [--ibl-sample] number of sample to compute a texel, default 1024"
        " -t [--thumbnail-size] width size of generated thumbnail, default 256, it will generate 256x128 thumbnail"
        " -d a debug flag to write intermediate texture"
        " -m [--no-threads] no multithreading, default is false";

    printf("%s", text);
}

int parseArgument(Options& options, int argc, char** argv)
{

    int opt;
    int optionIndex = 0;

    static const char* OPTSTR = "hdms:n:";
    static const struct option OPTIONS[] = {
        {"help", no_argument, nullptr, 'h'},           {"ibl-sample", required_argument, nullptr, 'n'},
        {"ibl-size", required_argument, nullptr, 's'}, {"no-threads", no_argument, nullptr, 's'},
        {"debug", no_argument, nullptr, 'd'},          {nullptr, 0, 0, 0} // termination of the option list
    };

    while ((opt = getopt_long(argc, argv, OPTSTR, OPTIONS, &optionIndex)) >= 0)
    {
        const char* arg = optarg ? optarg : "";
        if (opt == -1)
            break;

        switch (opt)
        {
        default:
        case 'h':
            printUsage();
            exit(0);
            break;
        case 'm':
            options.nbThreads = 1;
            break;
        case 's':
            options.cubemapSize = atoi(arg);
            break;
        case 'n':
            options.numSamples = atoi(arg);
            break;
        case 'd':
            options.debug = true;
            break;
        }
    }
    return optind;
}

int main(int argc, char** argv)
{
    Options options;
    int optionIndex = parseArgument(options, argc, argv);
    int numArgs = argc - optionIndex;
    if (numArgs < 1)
    {
        printUsage();
        return 1;
    }

    options.inputFilename = argv[optionIndex];
    if (numArgs > 1)
    {
        options.distDir = argv[optionIndex + 1];
    }

    printf("writing result to %s\n", options.distDir);
    if (makeDirectory(options.distDir))
    {
        printf("can't create director %s\n", options.distDir);
        return 1;
    }

    Image image;
    const char* distDir = options.distDir;

    int loadImageResult = envUtils::loadImage(image, options.inputFilename);
    if (loadImageResult == 0)
    {

        if (!options.nbThreads)
        {
            options.nbThreads = (int)std::thread::hardware_concurrency();
        }
        printf("using %d threads\n", options.nbThreads);

        // needs to be sure that the panorama is pow of 2
        Cubemap cm;
        envUtils::createCubemap(cm, options.cubemapSize);

        envUtils::clampImage(image, 255);

        envUtils::writeThumbnail(distDir, "thumbnail", image, options.thumbnailSize, options.thumbnailSize / 2);
        envUtils::equirectangularToCubemap(cm, image, options.nbThreads);

        if (options.debug)
            envUtils::writeCubemap_hdr(distDir, "input", cm);

        CubemapMipMap cmMipMap;
        envUtils::createCubemapMipMap(cmMipMap, cm);

        if (options.debug)
            envUtils::writeCubemapMipMapFaces_hdr(distDir, "mipmap", cmMipMap);

        CubemapMipMap cmPrefilter;
        envUtils::prefilterCubemapGGX(cmPrefilter, cmMipMap, options.numSamples, options.nbThreads);

        Image equirectangular;
        envUtils::packPrefilterCubemapToEquilateral(equirectangular, cmPrefilter, options.nbThreads);
        if (options.debug)
        {
            envUtils::writeImage_hdr(distDir, "equirectangular", equirectangular);
        }

        CubemapMipMap cmPrefilterFixUp;
        envUtils::resampleCubemap(cmPrefilterFixUp, cmPrefilter, options.nbThreads);

        if (options.debug)
            envUtils::writeCubemapMipMapFaces_hdr(distDir, "prefilter", cmPrefilter);

        Spherical spherical;
        envUtils::computeSphericalHarmonicsFromCubemap(spherical, cm);
        envUtils::writeSpherical_json(distDir, "spherical", spherical);
        envUtils::writeCubemapMipMap_luv(distDir, "prefilter", cmPrefilterFixUp);

        if (options.debug)
        {
            CubemapMipMap cmDecode;
            envUtils::readCubemapMipMap_luv(cmDecode, "test/prefilter.luv");
            envUtils::writeCubemapMipMapFaces_hdr(distDir, "prefilter-decode", cmDecode);
            envUtils::freeCubemapMipMap(cmDecode);
        }

        envUtils::writeCubemap_luv(distDir, "background", cmPrefilter.levels[2]);
        envUtils::writeImage_luv(distDir, "equirectangular", equirectangular);

        envUtils::freeCubemapMipMap(cmPrefilterFixUp);
        envUtils::freeCubemapMipMap(cmPrefilter);
        envUtils::freeCubemapMipMap(cmMipMap);
        envUtils::freeImage(equirectangular);
        envUtils::freeImage(image);
        envUtils::freeCubemap(cm);
    }

    return loadImageResult;
}
