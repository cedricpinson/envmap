#include "log.h"
#include "threadLines.h"

#include "Cubemap.h"
#include "CubemapMipMap.h"
#include "envUtils.h"
#include "utils.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

// define this to profile main functions that prefilter
// also needs to link with profile lib from gperftools
// then run pprof -top ./envmap ./prefilterCubemapGGX.txt
// ~/go/bin/pprof -callgrind ./prefilterCubemapGGX.txt  >callgrind.out
//#define PROFILER

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

//#define DEBUG_SAMPLE

namespace envUtils {

struct CacheEntry
{
    float direction[3];
    float lerp;
    unsigned char l0, l1;
    unsigned short padding;
};
struct CacheSample
{
    CacheEntry* samples;
    double totalWeight;
    int numSamples;
    int padding;

    void create(int num)
    {
        numSamples = num;
        samples = new CacheEntry[numSamples];
    }
    void free() { delete[] samples; }
};

// from filament I had the same but it's cleaner in filament
static void hemisphereImportanceSampleDggx(double3& H, double2 u, double a)
{
    const double phi = 2 * M_PI * u[0];
    // NOTE: (aa-1) == (a-1)(a+1) produces better fp accuracy
    const double cosTheta2 = (1 - u[1]) / (1 + (a + 1) * ((a - 1) * u[1]));
    const double cosTheta = sqrt(cosTheta2);
    const double sinTheta = sqrt(1 - cosTheta2);
    H[0] = sinTheta * cos(phi);
    H[1] = sinTheta * sin(phi);
    H[2] = cosTheta;
}

void precomputeSamples(CacheSample& cache, double roughnessLinear, size_t numSamples, size_t baseResolution)
{
    size_t tryNumSamples = numSamples;

    double3 H;
    // Solid angle covered by 1 pixel with 6 faces that are EnvMapSize X EnvMapSize
    const double omegaP = 4.0 * M_PI / (6.0 * baseResolution * baseResolution);

    const float maxLod = log2(baseResolution);

    // Original paper suggest biasing the mip to improve the results
    double mipBias = 1.0; // I tested that the result is better with bias 1

    double invSamples = 1.0 / numSamples;
    size_t count = 0;

    while (count != numSamples)
    {
        count = 0;

        // find the sequence to have desired sample hit NoL condition
        cache.totalWeight = 0.0;
        double hammersleyInvSamples = 1.0 / tryNumSamples;
        for (size_t s = 0; s < tryNumSamples && count < numSamples; s++)
        {

            // do the computation in local space and store the computed light vector L
            // in local space. It will be transformed in the main loop from tangent space and direction
            // https://placeholderart.wordpress.com/2015/07/28/implementation-notes-runtime-environment-map-filtering-for-image-based-lighting/

            /*
             *       (sampling)
             *            L         H (never calculated below)
             *            .\       /.
             *            . \     / .
             *            .  \   /  .
             *            .   \ /   .
             *         ---|----o----|-------> n
             *    cos(2*theta)    cos(theta)
             *       = n.L           = n.H
             *
             * Note: NoH == LoH
             * (H is the half-angle between L and V, and V == N)
             *
             * in local space: V == N == [0,0,1]
             */

            const double2 Xi = hammersley(s, hammersleyInvSamples);

            hemisphereImportanceSampleDggx(H, Xi, roughnessLinear);

#if 0
            // This produces the same result that the code below using the the non-simplified
            // equation. This let's us see that N == V and that L = -reflect(V, H)
            // Keep this for reference.
            const double3 N = {0, 0, 1};
            const double3 V = N;
            const double3 L = 2 * dot(H, V) * H - V;
            const double NoL = dot(N, L);
            const double NoH = dot(N, H);
            const double NoH2 = NoH * NoH;
            const double NoV = dot(N, V);
            const double LoH = dot(L, H);
#else
            // const double NoV = 1;
            const double NoH = H[2];
            const double NoH2 = NoH * NoH;
            const double NoL = 2 * NoH2 - 1;
            const double3 L(2 * NoH * H[0], 2 * NoH * H[1], NoL);
            // const double LoH = dot(L, H);
#endif

            if (NoL < 1.e-5)
                continue;

            CacheEntry& sample = cache.samples[count];

            // pre-filtered importance sampling
            // see: "Real-time Shading with Filtered Importance Sampling", Jaroslav Krivanek
            // see: "GPU-Based Importance Sampling, GPU Gems 3", Mark Colbert

            // Probability Distribution Function
            double Pdf = D_GGX(NoH, roughnessLinear) / 4.0;

            // Solid angle represented by this sample
            // float omegaS = 1.0 / (numSamples * Pdf);
            double omegaS = invSamples / Pdf;

            double mipLevel = fmin(fmax(0.5 * log2(omegaS / omegaP) + mipBias, 0.0), maxLod);

            sample.direction[0] = (float)L[0];
            sample.direction[1] = (float)L[1];
            sample.direction[2] = (float)L[2];

            double floorMipLevel = floor(mipLevel);
            sample.l0 = (unsigned char)floorMipLevel;
            sample.l1 = (unsigned char)ceil(mipLevel);
            sample.lerp = mipLevel - floorMipLevel;

            cache.totalWeight += NoL;

            count++;
        }

        tryNumSamples += numSamples - count;
    }
} // namespace envUtils

void writeWeightDistribution(const CacheSample& samples, double roughnessLinear, int numSamples, int mipLevel)
{
    Path path;
    char filename[32];
    snprintf(filename, 31, "sample_%d_%d.data", mipLevel, numSamples);
    createPath(path, "test", filename);
    FILE* fp = fopen(path, "w");
    fprintf(fp, "num samples %d mip level %d linear roughness %f\n", numSamples, mipLevel, roughnessLinear);
    double dir[3];
    for (int j = 0; j < numSamples; j++)
    {
        const CacheEntry& sample = samples.samples[j];
        // fprintf(fp, "%lf\n", samples.samples[j].direction[2]);
        float3ToDouble3(dir, sample.direction);
        fprintf(fp, "%f %f %f, %d %d %f\n", dir[0], dir[1], dir[2], sample.l0, sample.l1, (double)sample.lerp);
    }
    fclose(fp);
}

inline void computeBasisVector(float* tangentX, float* tangentY, const float* N)
{
    static float axisZ[3] = {0, 0, 1};
    static float axisX[3] = {1, 0, 0};
    const float* up = fabsf(N[2]) < 0.999f ? axisZ : axisX;

    float tmp[3];
    cross(tmp, up, N);
    normalize(tangentX, tmp);

    cross(tmp, N, tangentX);
    normalize(tangentY, tmp);
}

void getTrilinear(float* color, const Cubemap& cubemap0, const Cubemap& cubemap1, const float* direction,
                  float lerpFactor)
{
    Cubemap::Address address;
    Cubemap::getAddressFor(address, direction);

    float color0[3], color1[3];
    const Image& lod0 = cubemap0.faces[address.face];
    const Image& lod1 = cubemap1.faces[address.face];

    float x0 = address.s * lod0.width;
    float y0 = address.t * lod0.width;
    float x1 = address.s * lod1.width;
    float y1 = address.t * lod1.width;

    lod0.filterAt(color0, x0, y0);
    lod1.filterAt(color1, x1, y1);

    color[0] = color0[0] + (color1[0] - color0[0]) * lerpFactor;
    color[1] = color0[1] + (color1[1] - color0[1]) * lerpFactor;
    color[2] = color0[2] + (color1[2] - color0[2]) * lerpFactor;
}

struct PrefilterContext
{
    Cubemap* cubemapDest;
    const CubemapMipMap* cubemap;
    const CacheSample* samples;
    Cubemap::Face face;
    int padding;
    PrefilterContext(Cubemap* cm, Cubemap::Face face, const CubemapMipMap* cmMipMap, const CacheSample* sample)
        : cubemapDest(cm)
        , cubemap(cmMipMap)
        , samples(sample)
        , face(face)
    {}
};

// inline void rotateDirection(float3& L, float angle, const float3& l)
// {
//     float s, c, t;

//     s = sinf(angle);
//     c = cosf(angle);
//     t = 1.f - c;

//     L[0] = l[0] * c + l[1] * s;
//     L[1] = -l[0] * s + l[1] * c;
//     L[2] = l[2] * (t + c);
// }

inline void transformSampleLocal2World(float* world, const float* tangentX, const float* tangentY,
                                       const float* direction, const float* local)
{
    // lWorldSpace = tangentX * Lr[0] + tangentY * Lr[1] + direction * Lr[2]);
    world[0] = tangentX[0] * local[0] + tangentY[0] * local[1] + direction[0] * local[2];
    world[1] = tangentX[1] * local[0] + tangentY[1] * local[1] + direction[1] * local[2];
    world[2] = tangentX[2] * local[0] + tangentY[2] * local[1] + direction[2] * local[2];
}

void prefilterRangeLines(PrefilterContext context, int yStart, int yStop)
{
    const CacheSample& samples = *context.samples;
    Cubemap& cubemapDest = *context.cubemapDest;
    const CubemapMipMap& cubemap = *context.cubemap;
    Cubemap::Face faceIndex = context.face;

    int size = cubemapDest.size;
    Image& face = cubemapDest.faces[faceIndex];

    float tangentX[3];
    float tangentY[3];
    float direction[3];
    float color[3];
    float lWorldSpace[3];
    double prefilteredColor[3];

    for (int y = yStart; y <= yStop; ++y)
    {

        float* line = face.getPixel(0, y).ptr();
        for (int x = 0; x < size; x++)
        {
            prefilteredColor[0] = 0;
            prefilteredColor[1] = 0;
            prefilteredColor[2] = 0;

            cubemapDest.getDirectionFor(direction, faceIndex, x, y);

            computeBasisVector(tangentX, tangentY, direction);

            // https://placeholderart.wordpress.com/2015/07/28/implementation-notes-runtime-environment-map-filtering-for-image-based-lighting/
            // for the simplification

            // optimized lod version
            for (int i = 0; i < samples.numSamples; i++)
            {
                const CacheEntry& sample = samples.samples[i];
                const float* L = sample.direction;

                float NoL = L[2];

                transformSampleLocal2World(lWorldSpace, tangentX, tangentY, direction, L);

                const Cubemap& lod0 = cubemap.levels[sample.l0];
                const Cubemap& lod1 = cubemap.levels[sample.l1];

                getTrilinear(color, lod0, lod1, lWorldSpace, sample.lerp);

                prefilteredColor[0] += (double)(color[0] * NoL);
                prefilteredColor[1] += (double)(color[1] * NoL);
                prefilteredColor[2] += (double)(color[2] * NoL);
            }

            double invWeight = 1 / samples.totalWeight;

            line[x * 3 + 0] = (float)(prefilteredColor[0] * invWeight);
            line[x * 3 + 1] = (float)(prefilteredColor[1] * invWeight);
            line[x * 3 + 2] = (float)(prefilteredColor[2] * invWeight);
        }
    }
}
void prefilterCubemapGGX(CubemapMipMap& cmDst, const CubemapMipMap& cmSrc, size_t numSamples, int nbThread)
{
#ifdef PROFILER
    ProfilerStart("./prefilterCubemapGGX.txt");
#endif

    // build mipmap from cmSrc
    int size = cmSrc.levels[0].size;
    int maxMipMap = (int)log2(size);
    int numMipMap = maxMipMap + 1;

    // init dest result
    cmDst.init(numMipMap);
    for (int i = 0; i < numMipMap; i++)
    {
        envUtils::createCubemap(cmDst.levels[i], (int)pow(2, maxMipMap - i));
    }

    double stepRoughness = 1.0 / double(maxMipMap);

    char logString[24];
    snprintf(logString, 24, "prefiltering level %d", 0);
    auto t = logStart(logString);

    // copy for roughness = 0;
    memcpy(cmDst.levels[0].image.data, cmSrc.levels[0].image.data,
           cmSrc.levels[0].image.width * cmSrc.levels[0].image.height * sizeof(float3));
    cmDst.levels[0].makeSeamless();

    logEnd(t);

    CacheSample cacheSamples;
    cacheSamples.create(numSamples);

    for (int mipLevel = 1; mipLevel < numMipMap; mipLevel++)
    {
        snprintf(logString, 24, "prefiltering level %d", mipLevel);
        t = logStart(logString);
        double roughness = mipLevel * stepRoughness;
        double roughnessLinear = roughness * roughness;

        // frostbite, lagarde paper p67
        // http://www.frostbite.com/wp-content/uploads/2014/11/course_notes_moving_frostbite_to_pbr.pdf

        // precompute samples
        precomputeSamples(cacheSamples, roughnessLinear, numSamples, size);

#ifdef DEBUG_SAMPLE
        writeWeightDistribution(cacheSamples, roughnessLinear, numSamples, mipLevel);
#endif
        int nbLines = cmDst.levels[mipLevel].size;
        // iterate on each face to compute ggx
        for (int face = 0; face < 6; face++)
        {
            PrefilterContext context(&cmDst.levels[mipLevel], (Cubemap::Face)face, &cmSrc, &cacheSamples);
            threadLines(prefilterRangeLines, context, nbLines, nbThread);
        }
        logEnd(t);

        // needed for the remapping for webgl1 seamless cubemap
        cmDst.levels[mipLevel].makeSeamless();
    }

    cacheSamples.free();

#ifdef PROFILER
    ProfilerStop();
#endif
}

struct ResampleContext
{
    Cubemap* cubemapDest;
    const Cubemap* cubemap;
    Cubemap::Face face;
    int padding;
    ResampleContext(Cubemap* cm, Cubemap::Face face, const Cubemap* cmSrc)
        : cubemapDest(cm)
        , cubemap(cmSrc)
        , face(face)
    {}
};

// resample for seamless cubemap edge for webgl1 and more generally without extension seamless cubemap on gpu
// https://seblagarde.wordpress.com/2012/06/10/amd-cubemapgen-for-physically-based-rendering/
// http://code.google.com/p/nvidia-texture-tools/source/browse/trunk/src/nvtt/CubeSurface.cpp
void resampleCubemapRangeLines(ResampleContext context, int yStart, int yStop)
{
    Cubemap& cubemapDest = *context.cubemapDest;
    const Cubemap& cubemap = *context.cubemap;
    Cubemap::Face faceIndex = context.face;
    Cubemap::Address address;

    float size = cubemap.size;
    float direction[3];
    float color[3];
    float x0, y0;

    for (int y = yStart; y <= yStop; ++y)
    {

        float* line = cubemapDest.faces[faceIndex].getPixel(0, y).ptr();
        for (int x = 0; x < size; x++)
        {
            cubemap.getDirectionFixUpFor(direction, faceIndex, x, y);
            Cubemap::getAddressFor(address, direction);

            const Image& image = cubemap.faces[address.face];
            x0 = address.s * image.width;
            y0 = address.t * image.width;

            image.filterAt(color, x0, y0);

            *line++ = color[0];
            *line++ = color[1];
            *line++ = color[2];
        }
    }
}

void resampleCubemap(CubemapMipMap& cmDst, const CubemapMipMap& cmSrc, int nbThreads)
{
    auto t = logStart("resampling for seamless cubemap");

    // build mipmap from cmSrc
    int numMipMap = cmSrc.numLevel;
    cmDst.init(numMipMap);

    for (int i = 0; i < numMipMap; i++)
    {
        envUtils::createCubemap(cmDst.levels[i], cmSrc.levels[i].size);
    }

    for (int mipLevel = 0; mipLevel < numMipMap; mipLevel++)
    {
        const Cubemap* src = &cmSrc.levels[mipLevel];
        Cubemap* dst = &cmDst.levels[mipLevel];
        int nbLines = dst->size;

        // iterate on each face
        for (int face = 0; face < 6; face++)
        {
            ResampleContext context(dst, (Cubemap::Face)face, src);
            threadLines(resampleCubemapRangeLines, context, nbLines, nbThreads);
        }
    }

    logEnd(t);
}

} // namespace envUtils
