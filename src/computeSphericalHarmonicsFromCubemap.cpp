#include "Cubemap.h"
#include "Spherical.h"
#include "envUtils.h"
#include <stdio.h>
#include <string.h>

namespace envUtils {

// SH order use for approximation of irradiance cubemap is 5, mean 5*5 equals 25 coefficients
#define PI M_PI

// See Peter-Pike Sloan paper for these coefficients
static double SHBandFactor[NUM_SH_COEFFICIENT] = {1.0,         2.0 / 3.0,   2.0 / 3.0,   2.0 / 3.0,   1.0 / 4.0,
                                                  1.0 / 4.0,   1.0 / 4.0,   1.0 / 4.0,   1.0 / 4.0,   0.0,
                                                  0.0,         0.0,         0.0,         0.0,         0.0,
                                                  0.0, // The 4 band will be zeroed
                                                  -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0,
                                                  -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0};

void EvalSHBasis(const float4& dir, double* res)
{
    // Can be optimize by precomputing constant.
    static const double SqrtPi = sqrt(PI);

    double xx = (double)dir[0];
    double yy = (double)dir[1];
    double zz = (double)dir[2];

    // x[i] == pow(x, i), etc.
    double x[MAX_SH_ORDER + 1], y[MAX_SH_ORDER + 1], z[MAX_SH_ORDER + 1];
    x[0] = y[0] = z[0] = 1.;
    for (int i = 1; i < MAX_SH_ORDER + 1; ++i)
    {
        x[i] = xx * x[i - 1];
        y[i] = yy * y[i - 1];
        z[i] = zz * z[i - 1];
    }

    res[0] = (1 / (2. * SqrtPi));

    res[1] = -(sqrt(3 / PI) * yy) / 2.;
    res[2] = (sqrt(3 / PI) * zz) / 2.;
    res[3] = -(sqrt(3 / PI) * xx) / 2.;

    res[4] = (sqrt(15 / PI) * xx * yy) / 2.;
    res[5] = -(sqrt(15 / PI) * yy * zz) / 2.;
    res[6] = (sqrt(5 / PI) * (-1 + 3 * z[2])) / 4.;
    res[7] = -(sqrt(15 / PI) * xx * zz) / 2.;
    res[8] = sqrt(15 / PI) * (x[2] - y[2]) / 4.;

    res[9] = (sqrt(35 / (2. * PI)) * (-3 * x[2] * yy + y[3])) / 4.;
    res[10] = (sqrt(105 / PI) * xx * yy * zz) / 2.;
    res[11] = -(sqrt(21 / (2. * PI)) * yy * (-1 + 5 * z[2])) / 4.;
    res[12] = (sqrt(7 / PI) * zz * (-3 + 5 * z[2])) / 4.;
    res[13] = -(sqrt(21 / (2. * PI)) * xx * (-1 + 5 * z[2])) / 4.;
    res[14] = (sqrt(105 / PI) * (x[2] - y[2]) * zz) / 4.;
    res[15] = -(sqrt(35 / (2. * PI)) * (x[3] - 3 * xx * y[2])) / 4.;

    res[16] = (3 * sqrt(35 / PI) * xx * yy * (x[2] - y[2])) / 4.;
    res[17] = (-3 * sqrt(35 / (2. * PI)) * (3 * x[2] * yy - y[3]) * zz) / 4.;
    res[18] = (3 * sqrt(5 / PI) * xx * yy * (-1 + 7 * z[2])) / 4.;
    res[19] = (-3 * sqrt(5 / (2. * PI)) * yy * zz * (-3 + 7 * z[2])) / 4.;
    res[20] = (3 * (3 - 30 * z[2] + 35 * z[4])) / (16. * SqrtPi);
    res[21] = (-3 * sqrt(5 / (2. * PI)) * xx * zz * (-3 + 7 * z[2])) / 4.;
    res[22] = (3 * sqrt(5 / PI) * (x[2] - y[2]) * (-1 + 7 * z[2])) / 8.;
    res[23] = (-3 * sqrt(35 / (2. * PI)) * (x[3] - 3 * xx * y[2]) * zz) / 4.;
    res[24] = (3 * sqrt(35 / PI) * (x[4] - 6 * x[2] * y[2] + y[4])) / 16.;
}

// https://seblagarde.wordpress.com/2012/06/10/amd-cubemapgen-for-physically-based-rendering/
// void texelCoordToVectCubeMap(float* dirResult, int face, float ui, float vi, size_t size, int fixup)
// {

//     float u, v;

//     if (fixup)
//     {
//         // Code from Nvtt :
//         http://code.google.com/p/nvidia-texture-tools/source/browse/trunk/src/nvtt/CubeSurface.cpp

//         // transform from [0..res - 1] to [-1 .. 1], match up edges exactly.
//         u = (2.0f * ui / (size - 1.0f)) - 1.0f;
//         v = (2.0f * vi / (size - 1.0f)) - 1.0f;
//     }
//     else
//     {

//         // center ray on texel center
//         // generate a vector for each texel
//         u = (2.0f * (ui + 0.5f) / size) - 1.0f;
//         v = (2.0f * (vi + 0.5f) / size) - 1.0f;
//     }

//     float3 vecX = CubemapFace[face][0] * u;
//     float3 vecY = CubemapFace[face][1] * v;
//     float3 vecZ = CubemapFace[face][2];
//     float3 res = Vec3f(vecX + vecY + vecZ);
//     res.normalize();
//     dirResult[0] = res[0];
//     dirResult[1] = res[1];
//     dirResult[2] = res[2];
// }

// solid Angle
template <typename T> T AreaElement(const T x, const T y) { return atan2(x * y, sqrt(x * x + y * y + 1.0)); }

/** Original code from Ignacio Castaño
 * This formula is from Manne Öhrström's thesis.
 * Take two coordiantes in the range [-1, 1] that define a portion of a
 * cube face and return the area of the projection of that portion on the
 * surface of the sphere.
 **/
inline double texelPixelSolidAngle(const float& aU, const float& aV, const size_t width, const size_t height)
{
    // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
    // ( 0.5 is for texel center addressing)
    const double U = (2.0 * ((double)aU + 0.5) / width) - 1.0;
    const double V = (2.0 * ((double)aV + 0.5) / height) - 1.0;

    // Shift from a demi texel, mean 1.0 / size  with U and V in [-1..1]
    const double InvResolutionW = 1.0 / width;
    const double InvResolutionH = 1.0 / height;

    // U and V are the -1..1 texture coordinate on the current face.
    // Get projected area for this Texel
    const double x0 = U - InvResolutionW;
    const double y0 = V - InvResolutionH;
    const double x1 = U + InvResolutionW;
    const double y1 = V + InvResolutionH;
    const double SolidAngle = AreaElement(x0, y0) - AreaElement(x0, y1) - AreaElement(x1, y0) + AreaElement(x1, y1);

    return SolidAngle;
}

void computeSphericalHarmonicsFromCubemap(double* spherical, const Cubemap& cm)
{
    int size = cm.size;

    // First step - Generate SH coefficient for the diffuse convolution

    // allocate data for cubemap normalizer
    float4* cmNormalizeData = new float4[size * size * 6];
    float4* cmNormalize[6];
    for (int i = 0; i < 6; i++)
    {
        cmNormalize[i] = &cmNormalizeData[size * size * i];
    }

    float3 dir;

    // iterate over cube faces
    for (int iCubeFace = 0; iCubeFace < 6; iCubeFace++)
    {

        // fast texture walk, build normalizer cube map
        float4* texelPtr = cmNormalize[iCubeFace];

        for (int v = 0; v < size; v++)
        {
            for (int u = 0; u < size; u++)
            {
                // texelCoordToVectCubeMap(texelPtr, iCubeFace, (float)u, (float)v, size, 0);
                cm.getDirectionFor(dir, (Cubemap::Face)iCubeFace, u, v);
                double solidAngle = texelPixelSolidAngle(u, v, size, size);
                (*texelPtr)[0] = dir[0];
                (*texelPtr)[1] = dir[1];
                (*texelPtr)[2] = dir[2];
                (*texelPtr)[3] = solidAngle;
                texelPtr++;
            }
        }
    }

    // This is a custom implementation of D3DXSHProjectCubeMap to avoid to deal with LPDIRECT3DSURFACE9 pointer
    // Use Sh order 2 for a total of 9 coefficient as describe in http://www.cs.berkeley.edu/~ravir/papers/envmap/
    // accumulators are 64-bit floats in order to have the precision needed
    // over a summation of a large number of pixels
    double SHr[NUM_SH_COEFFICIENT];
    double SHg[NUM_SH_COEFFICIENT];
    double SHb[NUM_SH_COEFFICIENT];
    double SHdir[NUM_SH_COEFFICIENT];

    memset(SHr, 0, NUM_SH_COEFFICIENT * sizeof(double));
    memset(SHg, 0, NUM_SH_COEFFICIENT * sizeof(double));
    memset(SHb, 0, NUM_SH_COEFFICIENT * sizeof(double));
    memset(SHdir, 0, NUM_SH_COEFFICIENT * sizeof(double));

    double weightAccum = 0.0;

    for (int iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++)
    {
        for (int y = 0; y < size; y++)
        {
            float4* cmNormalizeFacePtr = &cmNormalize[iFaceIdx][(y * size)];
            const float3* cmFacePtr = &cm.faces[iFaceIdx].getPixel(0, y);

            for (int x = 0; x < size; x++)
            {
                // pointer to direction and solid angle in cube map associated with texel
                const float4& texelVect = cmNormalizeFacePtr[x];
                const double weight = (double)texelVect[3];

                // if (useSolidAngleWeighting)
                // { // solid angle stored in 4th channel of normalizer/solid angle cube map
                //     weight = texelVect[3];
                // }
                // else
                // { // all taps equally weighted
                //     weight = 1.0;
                // }

                EvalSHBasis(texelVect, SHdir);

                // Convert to double
                const float3& cmTexel = cmFacePtr[x];
                double R = (double)cmTexel[0];
                double G = (double)cmTexel[1];
                double B = (double)cmTexel[2];

                for (int i = 0; i < NUM_SH_COEFFICIENT; i++)
                {
                    SHr[i] += R * SHdir[i] * weight;
                    SHg[i] += G * SHdir[i] * weight;
                    SHb[i] += B * SHdir[i] * weight;
                }

                weightAccum += weight;
            }
        }
    }

    // Normalization - The sum of solid angle should be equal to the solid angle of the sphere (4 PI), so
    // normalize in order our weightAccum exactly match 4 PI.
    for (int i = 0; i < NUM_SH_COEFFICIENT; ++i)
    {
        SHr[i] *= 4.0 * PI / weightAccum;
        SHg[i] *= 4.0 * PI / weightAccum;
        SHb[i] *= 4.0 * PI / weightAccum;
    }

    for (int i = 0; i < NUM_SH_COEFFICIENT; ++i)
    {
        spherical[i * 3] = SHr[i] * SHBandFactor[i];
        spherical[i * 3 + 1] = SHg[i] * SHBandFactor[i];
        spherical[i * 3 + 2] = SHb[i] * SHBandFactor[i];
    }

#if 0
    // Second step - Generate cubemap from SH coefficient

    // regenerate normalization cubemap for the destination cubemap
    // clear pre-existing normalizer cube map
    // for(int iCubeFace=0; iCubeFace<6; iCubeFace++)
    // {
    //     normCubemap[iCubeFace].Clear();
    // }

    // Normalized vectors per cubeface and per-texel solid angle
    // BuildNormalizerSolidAngleCubemap(DstCubeImage->m_Width, m_NormCubeMap, a_FixupType);
    normCubemap.buildNormalizerSolidAngleCubemap(dstCubemap->getSize(), fixup);

    for (int iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++)
    {
        for (int y = 0; y < dstSize; y++)
        {
            normCubeRowStartPtr = &normCubemap.getImages().imageFace(iFaceIdx)[normCubeMapNumChannels * (y * dstSize)];
            dstCubeRowStartPtr = &dstCubemap->getImages().imageFace(iFaceIdx)[dstCubeMapNumChannels * (y * dstSize)];

            for (int x = 0; x < dstSize; x++)
            {
                // pointer to direction and solid angle in cube map associated with texel
                const float4& texelVect = &normCubeRowStartPtr[normCubeMapNumChannels * x];

                EvalSHBasis(texelVect, SHdir);

                // get color value
                float R = 0.0f, G = 0.0f, B = 0.0f;

                for (int i = 0; i < NUM_SH_COEFFICIENT; ++i)
                {
                    R += (float)(SHr[i] * SHdir[i] * SHBandFactor[i]);
                    G += (float)(SHg[i] * SHdir[i] * SHBandFactor[i]);
                    B += (float)(SHb[i] * SHdir[i] * SHBandFactor[i]);
                }

                dstCubeRowStartPtr[(dstCubeMapNumChannels * x) + 0] = R;
                dstCubeRowStartPtr[(dstCubeMapNumChannels * x) + 1] = G;
                dstCubeRowStartPtr[(dstCubeMapNumChannels * x) + 2] = B;
                if (dstCubeMapNumChannels > 3)
                {
                    dstCubeRowStartPtr[(dstCubeMapNumChannels * x) + 3] = 1.0f;
                }
            }
        }
    }
#endif

    delete[] cmNormalizeData;
}

} // namespace envUtils
