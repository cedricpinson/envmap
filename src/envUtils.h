#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "Image.h"

namespace envUtils {

// cubemap
void initCubemap(Cubemap& cm, const Image& image);
void createCubemap(Cubemap& cm, int size);
void equirectangularToCubemap(Cubemap& dest, const Image& equi);
void writeCubemap_hdr(const char* dir, const char* filename, const Cubemap& cm);
void freeCubemap(Cubemap& cm);
void downsampleCubemapLevelBoxFilter(Cubemap& dst, const Cubemap& src);

inline const char* getCubemapFaceName(Cubemap::Face face)
{
    switch (face)
    {
    case Cubemap::NX:
        return "face-NX";
        break;
    case Cubemap::PX:
        return "face-PX";
        break;
    case Cubemap::NY:
        return "face-NY";
        break;
    case Cubemap::PY:
        return "face-PY";
        break;
    case Cubemap::NZ:
        return "face-NZ";
        break;
    case Cubemap::PZ:
        return "face-PZ";
        break;
    }
    return "NONE";
}

// prefilterCubemapGGX
void prefilterCubemapGGX(CubemapMipMap& cmResult, const CubemapMipMap& cmSourceMipMap, size_t numSamples);

// cubemap mipmap
void createCubemapMipMap(CubemapMipMap& cmResultMipMap, const Cubemap& cmSource);
int writeCubemapMipMap_hdr(const char* dir, const char* basename, const CubemapMipMap& cm);
int writeCubemapMipMapFaces_hdr(const char* dir, const char* basename, const CubemapMipMap& cm);

// image
int loadImage(Image& image, const char* filename);
int writeImage_hdr(const char* filename, const Image& image);
void createImage(Image& image, int width, int height);
void freeImage(Image& image);
void clampImage(Image& src, float maxValue = 255);

// spherical harmonics coefficient
void computeSphericalHarmonicsFromCubemap(double* spherical, const Cubemap& cm);
int writeSpherical_json(const char* filename, double*);
} // namespace envUtils
