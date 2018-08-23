#include <atomic>
#include <iostream>
#include <thread>

#include "Cubemap.h"
#include "envUtils.h"
#include "utils.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

//#define DEBUG_SAMPLE

namespace envUtils {

struct CacheEntry
{
    float3 direction;
    float lerp;
    unsigned char l0, l1;
};
struct CacheSample
{
    CacheEntry* samples;
    double totalWeight;
    int numSamples;
    CacheSample(size_t numSamples)
        : numSamples(numSamples)
    {
        samples = new CacheEntry[numSamples];
    }
    ~CacheSample() { delete[] samples; }
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

void precomputeSamples(CacheSample& cache, float roughnessLinear, size_t numSamples, size_t baseResolution)
{
    size_t tryNumSamples = numSamples;

    double3 H;
    // Solid angle covered by 1 pixel with 6 faces that are EnvMapSize X EnvMapSize
    const double omegaP = 4.0 * M_PI / (6.0 * baseResolution * baseResolution);

    const float maxLod = log2(baseResolution);

    // Original paper suggest biasing the mip to improve the results
    float mipBias = 1.0f; // I tested that the result is better with bias 1

    float invSamples = 1.0 / numSamples;
    size_t count = 0;

    while (count != numSamples)
    {
        count = 0;

        // find the sequence to have desired sample hit NoL condition
        cache.totalWeight = 0.0;
        float hammersleyInvSamples = 1.0 / tryNumSamples;
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
            double Pdf = D_GGX(NoH, roughnessLinear) / 4.0f;

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

void writeWeightDistribution(const CacheSample& samples, float roughnessLinear, int numSamples, int mipLevel)
{
    Path path;
    char filename[32];
    snprintf(filename, 31, "sample_%d_%d.data", mipLevel, numSamples);
    create_path(path, "test", filename);
    FILE* fp = fopen(path, "w");
    fprintf(fp, "num samples %d mip level %d linear roughness %f\n", numSamples, mipLevel, roughnessLinear);
    for (int j = 0; j < numSamples; j++)
    {
        const CacheEntry& sample = samples.samples[j];
        // fprintf(fp, "%lf\n", samples.samples[j].direction[2]);
        fprintf(fp, "%f %f %f, %d %d %f\n", sample.direction[0], sample.direction[1], sample.direction[2], sample.l0,
                sample.l1, sample.lerp);
    }
    fclose(fp);
}

void computeBasisVector(float3& tangentX, float3& tangentY, const float3& N)
{
    static float3 axisZ = {0, 0, 1};
    static float3 axisX = {1, 0, 1};
    const float3& up = fabs(N[2]) < 0.999 ? axisZ : axisX;
    tangentX = normalize(cross(up, N));
    tangentY = normalize(cross(N, tangentX));
}

void getAddressFor(Cubemap::Address& addr, const float3& r)
{
    double sc, tc, ma;
    const float rx = fabs(r[0]);
    const float ry = fabs(r[1]);
    const float rz = fabs(r[2]);
    if (rx >= ry && rx >= rz)
    {
        ma = rx;
        if (r[0] >= 0)
        {
            addr.face = Cubemap::Face::PX;
            sc = -r[2];
            tc = -r[1];
        }
        else
        {
            addr.face = Cubemap::Face::NX;
            sc = r[2];
            tc = -r[1];
        }
    }
    else if (ry >= rx && ry >= rz)
    {
        ma = ry;
        if (r[1] >= 0)
        {
            addr.face = Cubemap::Face::PY;
            sc = r[0];
            tc = r[2];
        }
        else
        {
            addr.face = Cubemap::Face::NY;
            sc = r[0];
            tc = -r[2];
        }
    }
    else
    {
        ma = rz;
        if (r[2] >= 0)
        {
            addr.face = Cubemap::Face::PZ;
            sc = r[0];
            tc = -r[1];
        }
        else
        {
            addr.face = Cubemap::Face::NZ;
            sc = -r[0];
            tc = -r[1];
        }
    }
    // ma is guaranteed to be >= sc and tc
    addr.s = (sc / ma + 1) * 0.5f;
    addr.t = (tc / ma + 1) * 0.5f;
}

void getTrilinear(float3& color, const Cubemap cubemap0, const Cubemap cubemap1, const float3& direction,
                  float lerpFactor)
{
    Cubemap::Address address;
    getAddressFor(address, direction);

    float3 color0, color1;
    const Image& lod0 = cubemap0.faces[address.face];
    const Image& lod1 = cubemap1.faces[address.face];

    float x0, y0;
    float x1, y1;
    x0 = address.s * lod0.width;
    y0 = address.t * lod0.width;
    x1 = address.s * lod1.width;
    y1 = address.t * lod1.width;

    lod0.filterAt(color0, address.s * lod0.width, address.t * lod0.width);
    lod1.filterAt(color1, address.s * lod1.width, address.t * lod1.width);

    color = lerp(color0, color1, lerpFactor);
#ifdef DEBUG_SAMPLE
    printf("face %d l0 %f %f - l1 %f %f - color %f %f %f \n", address.face, x0, y0, x1, y1, color[0], color[1],
           color[2]);
#endif
}

void prefilterRangeLines(Cubemap* cubemapDestPtr, Cubemap::Face faceIndex, const CubemapMipMap* cubemapPtr,
                         const CacheSample* samplesPtr, int yStart, int yStop)
{
    const CacheSample& samples = *samplesPtr;
    Cubemap& cubemapDest = *cubemapDestPtr;
    const CubemapMipMap& cubemap = *cubemapPtr;

    int size = cubemapDest.size;
    Image& face = cubemapDest.faces[faceIndex];

    float3 tangentX;
    float3 tangentY;
    float3 direction;
    float3 color;
    float3 LworldSpace;

    for (int y = yStart; y <= yStop; ++y)
    {

        float3* line = &face.getPixel(0, y);
        for (int x = 0; x < size; x++)
        {

            cubemapDest.getDirectionFor(direction, faceIndex, x, y);

            double3 prefilteredColor = {0, 0, 0};

            computeBasisVector(tangentX, tangentY, direction);

            // see getPrecomputedLightInLocalSpace in Math
            // and
            // https://placeholderart.wordpress.com/2015/07/28/implementation-notes-runtime-environment-map-filtering-for-image-based-lighting/
            // for the simplification

            // optimized lod version
#ifdef DEBUG_SAMPLE
            printf("%f w %f face %d direction %f %f %f , %d %d\n", roughnessLinear, samples.totalWeight, faceIndex,
                   direction[0], direction[1], direction[2], x, y);
#endif
            for (int i = 0; i < samples.numSamples; i++)
            {
                const CacheEntry& sample = samples.samples[i];
                const float3& L = sample.direction;

                float NoL = L[2];
                float3 LworldSpace(tangentX * L[0] + tangentY * L[1] + direction * L[2]);

                const Cubemap& lod0 = cubemap.levels[sample.l0];
                const Cubemap& lod1 = cubemap.levels[sample.l1];

#ifdef DEBUG_SAMPLE
                printf("sample %d, %f %f %f (%d-%d) - ", i, LworldSpace[0], LworldSpace[1], LworldSpace[2], sample.l0,
                       sample.l1);
#endif
                getTrilinear(color, lod0, lod1, LworldSpace, sample.lerp);

                color *= NoL;

                prefilteredColor[0] += color[0];
                prefilteredColor[1] += color[1];
                prefilteredColor[2] += color[2];
            }

            prefilteredColor *= 1 / samples.totalWeight;
#ifdef DEBUG_SAMPLE
            printf("finale color %f %f %f\n", prefilteredColor[0], prefilteredColor[1], prefilteredColor[2]);
#endif

            line[x][0] = (float)prefilteredColor[0];
            line[x][1] = (float)prefilteredColor[1];
            line[x][2] = (float)prefilteredColor[2];
        }
    }
}
#define USE_THREAD

void prefilterFace(Cubemap& cubemapDest, Cubemap::Face faceIndex, const CubemapMipMap& cubemap, const CacheSample& samples, int nbThread)
{
    int nbLines = cubemapDest.size;
#ifdef USE_THREAD
    std::thread threadList[64];

    if (nbThread > 64)
        nbThread = 64;

    if (nbLines < nbThread)
    {
        nbThread = nbLines;
    }

    float linesPerThread = nbLines / nbThread;
    int startY = 0;
    int stopY;
    for (int i = 0; i < nbThread; i++)
    {
        stopY = startY + ceil(linesPerThread);
        if (stopY > nbLines - 1)
            stopY = nbLines - 1;

        threadList[i] = std::thread(prefilterRangeLines, &cubemapDest, faceIndex, &cubemap, &samples,
                                    startY, stopY);

        startY = stopY + 1;
    }

    for (int i = 0; i < nbThread; i++)
    {
        threadList[i].join();
    }
#else
    prefilterRangeLines(&cubemapDest, faceIndex, &cubemap, &samples, 0, nbLines - 1);
#endif
}

void prefilterCubemapGGX(CubemapMipMap& cmDst, const CubemapMipMap& cmSrc, size_t numSamples)
{
    int nbThread;
#ifdef USE_THREAD
    nbThread = std::thread::hardware_concurrency();
    printf("using %d threads\n", nbThread);
#endif

    // build mipmap from cmSrc
    int size = cmSrc.levels[0].size;
    int maxMipMap = log2(size);
    int numMipMap = maxMipMap + 1;

    // init dest result
    cmDst.init(numMipMap);
    for (int i = 0; i < numMipMap; i++)
    {
        envUtils::createCubemap(cmDst.levels[i], pow(2, maxMipMap - i));
    }

    float stepRoughness = 1.0 / float(maxMipMap);

    // copy for roughness = 0;
    memcpy(cmDst.levels[0].image.data, cmSrc.levels[0].image.data,
           cmSrc.levels[0].image.width * cmSrc.levels[0].image.height * sizeof(float3));

    CacheSample cacheSamples(numSamples);

    for (int mipLevel = 1; mipLevel < 2 // numMipMap
         ;
         mipLevel++)
    {
        float roughness = mipLevel * stepRoughness;
        float roughnessLinear = roughness * roughness;

        // frostbite, lagarde paper p67
        // http://www.frostbite.com/wp-content/uploads/2014/11/course_notes_moving_frostbite_to_pbr.pdf

        // precompute samples
        precomputeSamples(cacheSamples, roughnessLinear, numSamples, size);

        writeWeightDistribution(cacheSamples, roughnessLinear, numSamples, mipLevel);

        // iterate on each face to compute ggx
        for (int face = 0; face < 6; face++)
        {
            prefilterFace(cmDst.levels[mipLevel], (Cubemap::Face)face, cmSrc, cacheSamples, nbThread);
        }
    }
} // namespace envUtils

} // namespace envUtils
