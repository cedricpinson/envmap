#include "Cubemap.h"
#include "Image.h"

namespace envUtils {

    // cubemap
    void initCubemap(Cubemap& cm, const Image& image);
    void createCubemap(Cubemap& cm, int size);
    void equirectangularToCubemap(Cubemap& dest, const Image& equi);
    void writeCubemap_hdr(const char* dir, const Cubemap& cm);
    void freeCubemap(Cubemap& cm);

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

    // image
    int loadImage(Image& image, const char* filename);
    int writeImage_hdr(const char* filename, const Image& image);
    void createImage(Image& image, int width, int height);
    void freeImage(Image& image);

    // spherical harmonics coefficient
    void computeSphericalHarmonicsFromCubemap(double* spherical, const Cubemap& cm);
    int writeSpherical_json(const char* filename, double*);
}
