#pragma once

struct Cubemap;
struct CubemapMipMap;
struct Image;

namespace io {

int writeCubemapMipMapDebug(const char* dir, const char* basename, const CubemapMipMap& cm);
void writeCubemapDebug(const char* dir, const char* basename, const Cubemap& cm);

int writeCubemapMipMap_luv(const char* path, const CubemapMipMap& cm);
    //int writeCubemapMipMap_rgbm(const char* path, const CubemapMipMap& cm);
int readCubemapMipMap_luv(CubemapMipMap& cm, const char* path);

int writeCubemap_luv(const char* path, const Cubemap& cm);

int writeImage_hdr(const char* filename, const Image& image);
int writeImage_ldr(const char* filename, const Image& image);
int writeImage_luv(const char* filename, const Image& image);
    //int writeImage_rgbm(const char* filename, const Image& image);
int writeThumbnail(const char* dir, const char* basename, const Image& image, int width, int height);
int loadImage(Image& image, const char* filename);

int writeSpherical_json(const char* dir, const char* basename, double*);

} // namespace io
