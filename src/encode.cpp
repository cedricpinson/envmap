#include "encode.h"

// http://graphicrants.blogspot.fr/2009/04/rgbm-color-encoding.html
// M matrix, for encoding
const static float M[] = {0.2209, 0.3390, 0.4184, 0.1138, 0.6780, 0.7319, 0.0102, 0.1130, 0.2969};

// Inverse M matrix, for decoding
const static float InverseM[] = {6.0013, -2.700, -1.7995, -1.332, 3.1029, -5.7720, 0.3007, -1.088, 5.6268};

float3 mul(const float* v, const float* M)
{
    float3 result;

    result[0] = v[0] * M[0] + v[1] * M[3] + v[2] * M[6];
    result[1] = v[0] * M[1] + v[1] * M[4] + v[2] * M[7];
    result[2] = v[0] * M[2] + v[1] * M[5] + v[2] * M[8];

    return result;
}

float4 LogLuvEncode(const float* vRGB)
{
    float4 vResult;
    float3 Xp_Y_XYZp;
    float3 test;

    Xp_Y_XYZp = mul(vRGB, M);

    Xp_Y_XYZp = Xp_Y_XYZp.max(1e-6f);

    vResult[0] = Xp_Y_XYZp[0] / Xp_Y_XYZp[2];
    vResult[1] = Xp_Y_XYZp[1] / Xp_Y_XYZp[2];

    float Le = 2.0 * log2(Xp_Y_XYZp[1]) + 127.0;
    vResult[3] = frac(Le);
    vResult[2] = (Le - (floor(vResult[3] * 255.0f)) / 255.0f) / 255.0f;

    return vResult;
}

void encodeLUV(uint8_t* luvDst, const float* rgbSrc)
{
    float4 result = LogLuvEncode(rgbSrc);
    luvDst[0] = uint8_t(result[0] * 255.0);
    luvDst[1] = uint8_t(result[1] * 255.0);
    luvDst[2] = uint8_t(result[2] * 255.0);
    luvDst[3] = uint8_t(result[3] * 255.0);
}

float uncharted2Tonemap(float x)
{
    float A = 0.15;
    float B = 0.50;
    float C = 0.10;
    float D = 0.20;
    float E = 0.02;
    float F = 0.30;

    return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

void tonemap(uint8_t* rgbDst, const float* rgbSrc)
{
    float W = 11.2;

    float rgb[3];

    float exposureBias = 2.0;
    rgb[0] = uncharted2Tonemap(rgbSrc[0] * exposureBias);
    rgb[1] = uncharted2Tonemap(rgbSrc[1] * exposureBias);
    rgb[2] = uncharted2Tonemap(rgbSrc[2] * exposureBias);

    float whiteScale = 1.0 / uncharted2Tonemap(W);

    rgb[0] *= whiteScale;
    rgb[1] *= whiteScale;
    rgb[2] *= whiteScale;

    float gammaInv = 1.0 / 2.2;
    rgbDst[0] = (uint8_t)floor(255 * clamp(powf(rgb[0], gammaInv), 0.0f, 1.0f));
    rgbDst[1] = (uint8_t)floor(255 * clamp(powf(rgb[1], gammaInv), 0.0f, 1.0f));
    rgbDst[2] = (uint8_t)floor(255 * clamp(powf(rgb[2], gammaInv), 0.0f, 1.0f));
}

float3 LogLuvDecode(const float* vLogLuv)
{
    float Le = vLogLuv[2] * 255.0 + vLogLuv[3];
    float3 Xp_Y_XYZp;
    Xp_Y_XYZp[1] = exp2((Le - 127.0) / 2.0);
    Xp_Y_XYZp[2] = Xp_Y_XYZp[1] / vLogLuv[1];
    Xp_Y_XYZp[0] = vLogLuv[0] * Xp_Y_XYZp[2];
    float3 vRGB;

    vRGB = mul(&Xp_Y_XYZp[0], InverseM);

    return vRGB.max(0.0);
}

void decodeLUV(float* rgb, const uint8_t* luv)
{
    float4 value;
    value[0] = luv[0] * 1.0 / 255.0;
    value[1] = luv[1] * 1.0 / 255.0;
    value[2] = luv[2] * 1.0 / 255.0;
    value[3] = luv[3] * 1.0 / 255.0;

    float3 result = LogLuvDecode(value.ptr());

    rgb[0] = result[0];
    rgb[1] = result[1];
    rgb[2] = result[2];
}

float RGBMMaxRange = 8.0;

// https://gist.github.com/aras-p/1199797
void encodeRGBM(float rgb[3], uint8_t rgbm[4])
{

    // in our case,
    const float kRGBMMaxRange = RGBMMaxRange;
    const float kOneOverRGBMMaxRange = 1.0f / kRGBMMaxRange;

    // encode to RGBM, c = ARGB colors in 0..1 floats

    float r = rgb[0] * kOneOverRGBMMaxRange;
    float g = rgb[1] * kOneOverRGBMMaxRange;
    float b = rgb[2] * kOneOverRGBMMaxRange;

    float a = fmax(fmax(r, g), fmax(b, 1e-6f));
    a = ceilf(a * 255.0f) / 255.0f;

    rgbm[0] = uint8_t(fmin(r / a, 1.0f) * 255);
    rgbm[1] = uint8_t(fmin(g / a, 1.0f) * 255);
    rgbm[2] = uint8_t(fmin(b / a, 1.0f) * 255);
    rgbm[3] = uint8_t(fmin(a, 1.0f) * 255.0);
}

void decodeRGBM(uint8_t rgbm[4], float rgb[3])
{
    rgb[0] = rgbm[0] / 255.0 * RGBMMaxRange * rgbm[3] / 255.0;
    rgb[1] = rgbm[1] / 255.0 * RGBMMaxRange * rgbm[3] / 255.0;
    rgb[2] = rgbm[2] / 255.0 * RGBMMaxRange * rgbm[3] / 255.0;
}
