#include "Cubemap.h"
#include "envUtils.h"
#include "log.h"
#include "threadLines.h"

namespace envUtils {

struct EquirectangularProcessContext
{
    Cubemap* dst;
    const Image* src;
    Cubemap::Face face;
    int padding;
    EquirectangularProcessContext(Cubemap* dst, Cubemap::Face face, const Image* src)
        : dst(dst)
        , src(src)
        , face(face)
    {}
};

// this part is from filament https://github.com/google/filament/blob/master/tools/cmgen/src/CubemapUtils.cpp
void equirectangularToCubemapLines(EquirectangularProcessContext context, int startY, int stopY)
{
    Cubemap& dst = *context.dst;
    Cubemap::Face faceIdex = context.face;
    Image& face = dst.faces[faceIdex];
    const Image& src = *context.src;

    const int width = src.width;
    const int height = src.height;
    const double r = width * 0.5 * M_1_PI;
    int dim = face.width;

    for (int y = startY; y <= stopY; y++)
    {
        float3* data = &face.getPixel(0, y);
        for (int x = 0; x < dim; ++x, ++data)
        {
            // calculate how many samples we need based on dx, dy in the source
            // x =  cos(phi) sin(theta)
            // y = -sin(phi)
            // z = -cos(phi) cos(theta)
            float3 s0f;
            dst.getDirectionFor(s0f.ptr(), faceIdex, x, y);
            double3 s0 = s0f.toDouble();

            const double t0 = atan2(s0[0], -s0[2]);
            const double p0 = asin(s0[1]);
            float3 s1f;
            dst.getDirectionFor(s1f.ptr(), faceIdex, x + 1, y + 1);
            double3 s1 = s1f.toDouble();

            const double t1 = atan2(s1[0], -s1[2]);
            const double p1 = asin(s1[1]);
            const double dt = abs(t1 - t0);
            const double dp = abs(p1 - p0);
            const double dx = abs(r * dt);
            const double dy = abs(r * dp * s0[1]);
            const size_t numSamples = (size_t const)ceil(fmax(dx, dy));
            const double iNumSamples = 1.0 / numSamples;

            double3 c = double3(0, 0, 0);
            for (size_t sample = 0; sample < numSamples; sample++)
            {
                // Generate numSamples in our destination pixels and map them to input pixels
                const double2 h = hammersley(size_t(sample), iNumSamples);
                float3 sf;
                dst.getDirectionFor(sf.ptr(), faceIdex, float(x + h[0]), float(y + h[1]));
                double3 s = sf.toDouble();
                double xf = atan2(s[0], -s[2]) * M_1_PI; // range [-1.0, 1.0]
                double yf = asin(-s[1]) * (2 * M_1_PI);  // range [-1.0, 1.0]
                xf = (xf + 1) * 0.5 * (width - 1);       // range [0, width [
                yf = (yf + 1) * 0.5 * (height - 1);      // range [0, height[
                // we can't use filterAt() here because it reads past the width/height
                // which is okay for cubmaps but not for square images

                int xSample = (int)xf;
                int ySample = (int)yf;
                xSample = xSample < width ? xSample : width - 1;
                ySample = ySample < height ? ySample : height - 1;
                const float3& pixel = src.getPixel(xSample, ySample);
                c[0] += (double)pixel[0];
                c[1] += (double)pixel[1];
                c[2] += (double)pixel[2];
            }
            c *= iNumSamples;
            float3& resultPixel = *data;
            resultPixel[0] = c[0];
            resultPixel[1] = c[1];
            resultPixel[2] = c[2];
        }
    }
}

void equirectangularToCubemap(Cubemap& dst, const Image& src, int nbThread)
{
    auto t = logStart("equirectangularToCubemap");
    for (int f = 0; f < 6; f++)
    {
        EquirectangularProcessContext context(&dst, (Cubemap::Face)f, &src);
        threadLines(equirectangularToCubemapLines, context, dst.size, nbThread);
    }
    logEnd(t);
}

struct CubemapToEquirectangularProcessContext
{
    Image* dst;
    const Cubemap* src;
    CubemapToEquirectangularProcessContext(Image* dst, const Cubemap* src)
        : dst(dst)
        , src(src)
    {}
};

void cubemapToEquirectangularLines(CubemapToEquirectangularProcessContext context, int startY, int stopY)
{
    Image& dst = *context.dst;
    const Cubemap& cubemap = *context.src;

    int width = dst.width;
    int height = dst.height;

    float direction[3];
    float u, v, theta, phi;

    float PI_2 = 2.0 * M_PI;
    float invWidth = 1.0 / width;
    float invHeight = 1.0 / height;
    float sinTheta;

    Cubemap::Address address;
    float x0, y0;

    for (int y = startY; y <= stopY; y++)
    {
        v = (y + 0.5f) * invHeight;
        theta = v * (float)M_PI;

        float* data = dst.getPixel(0, y).ptr();
        for (int x = 0; x < width; ++x, data += 3)
        {
            u = (x + 0.5f) * invWidth;
            phi = u * PI_2;

            sinTheta = sinf(theta);
            direction[0] = -sinf(phi) * sinTheta;
            direction[1] = cos(theta);
            direction[2] = -cosf(phi) * sinTheta;

            Cubemap::getAddressFor(address, direction);

            const Image& image = cubemap.faces[address.face];
            x0 = address.s * cubemap.size;
            y0 = address.t * cubemap.size;
            image.filterAt(data, x0, y0);
        }
    }
}

void cubemapToEquirectangular(Image& dst, const Cubemap& src, int nbThread)
{
    for (int f = 0; f < 6; f++)
    {
        CubemapToEquirectangularProcessContext context(&dst, &src);
        threadLines(cubemapToEquirectangularLines, context, dst.height, nbThread);
    }
}

} // namespace envUtils
