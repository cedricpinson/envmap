#pragma once

#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "Image.h"

namespace envUtils {

// cubemap
void initCubemap(Cubemap& cm, const Image& image);
void createCubemap(Cubemap& cm, int size);
void equirectangularToCubemap(Cubemap& dest, const Image& equi, int nbThread);
void writeCubemap_hdr(const char* dir, const char* filename, const Cubemap& cm);
int writeCubemap_luv(const char* dir, const char* basename, const Cubemap& cm);
void freeCubemap(Cubemap& cm);
void downsampleCubemapLevelBoxFilter(Cubemap& dst, const Cubemap& src);

// prefilterCubemapGGX
void prefilterCubemapGGX(CubemapMipMap& cmResult, const CubemapMipMap& cmSourceMipMap, size_t numSamples, int nbThread);

// cubemap mipmap
void createCubemapMipMap(CubemapMipMap& cmResultMipMap, const Cubemap& cmSource);
int writeCubemapMipMap_hdr(const char* dir, const char* basename, const CubemapMipMap& cm);
int writeCubemapMipMapFaces_hdr(const char* dir, const char* basename, const CubemapMipMap& cm);
int writeCubemapMipMap_luv(const char* dir, const char* basename, const CubemapMipMap& cm);
int readCubemapMipMap_luv(CubemapMipMap& cm, const char* path);

// image
int loadImage(Image& image, const char* filename);
int writeImage_hdr(const char* filename, const Image& image);
void createImage(Image& image, int width, int height);
void freeImage(Image& image);
void clampImage(Image& src, float maxValue = 255);
int writeThumbnail(const char* dir, const char* basename, const Image& image, int width, int height);

// spherical harmonics coefficient
void computeSphericalHarmonicsFromCubemap(double* spherical, const Cubemap& cm);
int writeSpherical_json(const char* dir, const char* basename, double*);

} // namespace envUtils
