#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "Image.h"
#include "Light.h"
#include "Spherical.h"
#include "envUtils.h"
#include "io.h"
#include "lightExtraction.h"
#include "package.h"
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
        "envmap is a tool to precompute environment for PBR rendering\n"
        "usages:\n"
        " envmap [options] panorama_env.hdr\n"
        "\n"
        "example\n"
        " envmap -s 512 -n 4096 panorama_env.hdr\n"
        "\n"
        "Options:\n"
        " -s [--ibl-size] to define the cubemap resolution computed, default 256\n"
        " -n [--ibl-sample] number of sample to compute a texel, default 1024\n"
        " -t [--thumbnail-size] width size of generated thumbnail, default 256, it will generate 256x128 thumbnail\n"
        " -d a debug flag to write intermediate texture\n"
        " -m [--no-threads] no multithreading, default is false\n";

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
    Path path;

    int loadImageResult = io::loadImage(image, options.inputFilename);
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

        Image thumbnail;
        envUtils::resizeImage(thumbnail, image, options.thumbnailSize, options.thumbnailSize / 2);

        envUtils::equirectangularToCubemap(cm, image, options.nbThreads);

        if (options.debug)
        {
            io::writeCubemapDebug(distDir, "input", cm);
        }

        CubemapMipMap cmMipMap;
        envUtils::createCubemapMipMap(cmMipMap, cm);

        if (options.debug)
            io::writeCubemapMipMapDebug(distDir, "mipmap", cmMipMap);

        CubemapMipMap cmPrefilter;
        envUtils::prefilterCubemapGGX(cmPrefilter, cmMipMap, options.numSamples, options.nbThreads);

        Image equirectangular;
        envUtils::packPrefilterCubemapToEquilateral(equirectangular, cmPrefilter, options.nbThreads);
        if (options.debug)
        {
            createPath(path, distDir, "equirectangular.hdr");
            io::writeImage_hdr(path, equirectangular);
        }

        CubemapMipMap cmPrefilterFixUp;
        envUtils::resampleCubemap(cmPrefilterFixUp, cmPrefilter, options.nbThreads);

        if (options.debug)
        {
            io::writeCubemapMipMapDebug(distDir, "prefilter", cmPrefilter);
        }

        Spherical spherical;
        envUtils::computeSphericalHarmonicsFromCubemap(spherical, cm);

        if (options.debug)
        {
            CubemapMipMap cmDecode;
            io::readCubemapMipMap_luv(cmDecode, "test/prefilter.luv");
            io::writeCubemapMipMapDebug(distDir, "prefilter-decode", cmDecode);
            envUtils::freeCubemapMipMap(cmDecode);
        }

        Light light;
        Image inputLightExtraction;
        // io::loadImage(inputLightExtraction, "extraction_input.exr");
        envUtils::resizeImage(inputLightExtraction, image, 1024, 512);
        io::writeImage_hdr("extraction_input-0.hdr", inputLightExtraction);
        if (extractMainLight(light, inputLightExtraction) != 0)
        {
            printf("Did not find a main light from the environment\n");
        }
        envUtils::freeImage(inputLightExtraction);

        pkg::Package package(distDir);

        package.setSpherical(&spherical);
        package.addThumbnail(thumbnail);
        package.addPrefilterCubemap(cmPrefilterFixUp, pkg::ImageEncoding::luv);
        package.addPrefilterEquirectangular(equirectangular, pkg::ImageEncoding::luv);
        package.addBackground(cmPrefilter.levels[2], pkg::ImageEncoding::luv);
        package.addBackground(cmPrefilter.levels[0], pkg::ImageEncoding::luv);
        package.addBackground(cmPrefilter.levels[4], pkg::ImageEncoding::luv);

        package.addPrefilterCubemap(cmPrefilterFixUp, pkg::ImageEncoding::rgbm);
        package.addPrefilterEquirectangular(equirectangular, pkg::ImageEncoding::rgbm);
        package.addBackground(cmPrefilter.levels[2], pkg::ImageEncoding::rgbm);
        package.addBackground(cmPrefilter.levels[0], pkg::ImageEncoding::rgbm);
        package.addBackground(cmPrefilter.levels[4], pkg::ImageEncoding::rgbm);

        package.write();

        envUtils::freeCubemapMipMap(cmPrefilterFixUp);
        envUtils::freeCubemapMipMap(cmPrefilter);
        envUtils::freeCubemapMipMap(cmMipMap);
        envUtils::freeImage(equirectangular);
        envUtils::freeImage(thumbnail);
        envUtils::freeImage(image);
        envUtils::freeCubemap(cm);
    }

    return loadImageResult;
}
