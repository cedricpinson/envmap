#include "Cubemap.h"
#include "envUtils.h"
#include "utils.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

namespace envUtils {

struct CacheSample
{
    float4* samples;
    double totalWeight;
    CacheSample(size_t maxSample) { samples = new float4[maxSample]; }
    ~CacheSample() { delete[] samples; }
};

bool computeLightSampleInLocalSpace(float4& result, size_t i, float invNumSamples, size_t size, float roughnessLinear)
{
    // do the computation in local space and store the computed light vector L
    // in local space. It will be transformed in the main loop from tangent space

    // https://placeholderart.wordpress.com/2015/07/28/implementation-notes-runtime-environment-map-filtering-for-image-based-lighting/
    // in local space view == normal == 0,0,1
    static const float3 V = {0, 0, 1};

    float2 Xi = hammersley(i, invNumSamples);

    float roughness = roughnessLinear * roughnessLinear;
    float Phi = 2.0 * M_PI * Xi[0];
    float CosTheta = sqrt((1.0 - Xi[1]) / (1.0 + (roughness * roughness - 1.0) * Xi[1]));
    float SinTheta = sqrt(1.0 - CosTheta * CosTheta);
    float3 H;
    H[0] = SinTheta * cos(Phi);
    H[1] = SinTheta * sin(Phi);
    H[2] = CosTheta;
    H.normalize();

    float3 L = H * (dot(V, H) * 2.0) - V;
    L.normalize();

    float NoL = L[2];

    if (NoL <= 0.0f)
        return false;

    // H[2] is NoH in local space
    // adds 1.e-5 to avoid D_GGX / 0.0
    float NoH = H[2] + 1.e-5;

    // data to evaluate pdf
    float VoH = NoH;

    // Probability Distribution Function
    float Pdf = D_GGX(NoH, roughnessLinear) * NoH / (4.0f * VoH);

    // Solid angle represented by this sample
    // float omegaS = 1.0 / (numSamples * Pdf);
    float omegaS = invNumSamples / Pdf;

    // Solid angle covered by 1 pixel with 6 faces that are EnvMapSize X EnvMapSize
    const float omegaP = 4.0 * M_PI / (6.0 * size * size);

    // Original paper suggest biasing the mip to improve the results
    float mipBias = 1.0f; // I tested that the result is better with bias 1
    double maxLod = log2(size);
    float mipLevel = fmin(fmax(0.5 * log2(omegaS / omegaP) + mipBias, 0.0), maxLod);

    result[0] = L[0];
    result[1] = L[1];
    result[2] = L[2];
    result[3] = mipLevel;
    return true;
}

void precomputeSamples(CacheSample& cache, float roughnessLinear, size_t numSamples, size_t baseResolution)
{
    size_t tryNumSamples = numSamples;
    size_t count = 0;
    while (count < numSamples)
    {
        count = 0;
        size_t index = 0;

        // find the sequence to have desired sample hit NoL condition
        cache.totalWeight = 0.0;
        float invNumSamples = 1.0 / tryNumSamples;
        for (size_t s = 0; s < tryNumSamples; s++)
        {

            if (computeLightSampleInLocalSpace(cache.samples[index], s, invNumSamples, baseResolution, roughnessLinear))
            {
                cache.totalWeight += cache.samples[index][2];
                count++;
                index++;
            }
        }

        tryNumSamples += numSamples - count;
    }
}

void prefilterCubemapGGX(CubemapMipMap& cmDst, const CubemapMipMap& cmSrc, size_t numSamples)
{
    // build mipmap from cmSrc
    size_t size = cmSrc.levels[0].size;
    size_t numMipMap = log2(size);

    // init dest result
    cmDst.init(numMipMap);
    for (int i = 0; i < numMipMap; i++)
    {
        envUtils::createCubemap(cmDst.levels[i], pow(2, numMipMap - i));
    }

    float stepRoughness = 1.0 / float(numMipMap - 1);

    // copy for roughness = 0;
    memcpy(cmDst.levels[0].image.data, cmSrc.levels[0].image.data,
           cmSrc.levels[0].image.width * cmSrc.levels[0].image.height * sizeof(float3));

    CacheSample cacheSamples(numSamples);

    for (int i = 1; i < numMipMap; i++)
    {
        float roughness = i * stepRoughness;
        float roughnessLinear = roughness * roughness;

        // frostbite, lagarde paper p67
        // http://www.frostbite.com/wp-content/uploads/2014/11/course_notes_moving_frostbite_to_pbr.pdf

        // precompute samples
        precomputeSamples(cacheSamples, roughnessLinear, numSamples, size);
        // iterate on each face to compute ggx

        Path path;
        char filename[32];
        snprintf(filename, 31, "sample_%d_%d.data", i, (int)numSamples);
        create_path(path, "test", filename);
        FILE* fp = fopen(path, "w");
        fprintf(fp, "num samples %d mip level %d linear roughness %f\n", (int)numSamples, i, roughnessLinear);
        for (int j = 0; j < numSamples; j++)
        {
            fprintf(fp, "%f\n", cacheSamples.samples[j][2]);
        }
        fclose(fp);
    }
}

} // namespace envUtils
