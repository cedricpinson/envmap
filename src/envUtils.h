#pragma once

#include <sys/types.h>

struct CubemapMipMap;
struct Cubemap;
struct Image;

namespace envUtils {

// cubemap
void initCubemap(Cubemap& cm, const Image& image);
void createCubemap(Cubemap& cm, int size);
void equirectangularToCubemap(Cubemap& dest, const Image& equi, int nbThread);
void cubemapToEquirectangular(Image& dst, const Cubemap& src, int nbThread);
void freeCubemap(Cubemap& cm);
void downsampleCubemapLevelBoxFilter(Cubemap& dst, const Cubemap& src);

// prefilterCubemapGGX
void prefilterCubemapGGX(CubemapMipMap& cmResult, const CubemapMipMap& cmSourceMipMap, size_t numSamples, int nbThread);
void resampleCubemap(CubemapMipMap& cmDst, const CubemapMipMap& cmSrc, int nbThreads);
void packPrefilterCubemapToEquilateral(Image& equirectangular, const CubemapMipMap& src, int nbThreads);

// cubemap mipmap
void createCubemapMipMap(CubemapMipMap& cmResultMipMap, const Cubemap& cmSource);
void freeCubemapMipMap(CubemapMipMap& cmMipMap);

// image
void createImage(Image& image, int width, int height);
void freeImage(Image& image);
void clampImage(Image& src, float maxValue = 255);
void resizeImage(Image& dst, const Image& image, int width, int height);

// spherical harmonics coefficient
void computeSphericalHarmonicsFromCubemap(double* spherical, const Cubemap& cm);

} // namespace envUtils
