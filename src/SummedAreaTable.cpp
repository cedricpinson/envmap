#include "SummedAreaTable.h"
#include "Image.h"
#include <assert.h>

void SummedAreaTable::free() { delete[] _sat; }
void SummedAreaTable::create(const Image& image)
{
    int width = image.width;
    int height = image.height;
    _width = width;
    _height = height;

    const int imgSize = width * height;

    _sat = new double[imgSize];

    double weightAccum = 0.0;

    // solid angle for 1 pixel on equi map
    double weight = (4.0 * M_PI) / ((double)(imgSize));

    double minLum = DBL_MAX;
    double maxLum = DBL_MIN;
    double minPonderedLum = DBL_MAX;
    double maxPonderedLum = DBL_MIN;
    _sumLuminance = 0.0;

    for (int y = 0; y < height; ++y)
    {

        const double posY = (double)(y + 1.0) / (double)(height + 1.0);

        // the latitude-longitude format overrepresents the area of regions near the poles.
        // To compensate for this, the pixels of the probe image
        // should first be scaled by cosφ.
        // (φ == 0 at middle height of image input)
        const double solidAngle = cos(M_PI * (posY - 0.5)) * weight;

        for (int x = 0; x < width; ++x)
        {
            const int i = y * width + x;

            double3 pixel;
            float3ToDouble3(pixel.ptr(), image.getPixel(x, y).ptr());

            double ixy = luminance(pixel);

            // update Min/Max before pondering
            minLum = fmin(ixy, minLum);
            maxLum = fmax(ixy, maxLum);

            pixel *= solidAngle * imgSize;

            // pondering luminance for unpondered colors makes more sense
            ixy *= solidAngle;

            _sat[i] = ixy;

            weightAccum += solidAngle;
            _sumLuminance += ixy;
        }
    }
    // store for later use.
    _weightAccum = weightAccum;

    // normalize in order our image Accumulation exactly match 4 PI.
    const double normalizer = (4.0 * M_PI) / weightAccum;

    _sumLuminance *= normalizer;

    for (int i = 0; i < imgSize; ++i)
    {
        _sat[i] *= normalizer;
        minPonderedLum = fmin(_sat[i], minPonderedLum);
        maxPonderedLum = fmax(_sat[i], maxPonderedLum);
    }

    // enhances precision of SAT
    // make values be around [0.0, 0.5]
    // https://developer.amd.com/wordpress/media/2012/10/SATsketch-siggraph05.pdf
    const double rangePonderedLum = maxPonderedLum - minPonderedLum;

    for (int i = 0; i < imgSize; ++i)
    {
        _sat[i] = ((_sat[i] - minPonderedLum) / rangePonderedLum) * 0.5;
    }

    // now we sum
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            const int i = y * width + x;
            // https://en.wikipedia.org/wiki/Summed_area_table
            _sat[i] = _sat[i] + I(x - 1, y) + I(x, y - 1) - I(x - 1, y - 1);
        }
    }
}
