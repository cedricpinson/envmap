#include "lightExtraction.h"
#include "Image.h"
#include "Light.h"
#include "SummedAreaTable.h"
#include <stdio.h>
#include <stdlib.h>

int lightInvCmpSum(const void* a0, const void* b0)
{
    const Light& a = *static_cast<const Light*>(a0);
    const Light& b = *static_cast<const Light*>(b0);
    double v = b._sum - a._sum;
    if (v > 0.0)
        return 1;
    else if (v < 0.0)
        return -1;
    return 0;
}

int lightCmpAreaSize(const void* a0, const void* b0)
{
    const Light& a = *static_cast<const Light*>(a0);
    const Light& b = *static_cast<const Light*>(b0);
    double v = a._areaSize - b._areaSize;
    if (v > 0.0)
        return 1;
    else if (v < 0.0)
        return -1;
    return 0;
}

int lightCmpLumAverage(const void* a0, const void* b0)
{
    const Light& a = *static_cast<const Light*>(a0);
    const Light& b = *static_cast<const Light*>(b0);
    double v = a._lumAverage - b._lumAverage;
    if (v > 0.0)
        return 1;
    else if (v < 0.0)
        return -1;
    return 0;
}

/**
 * A subregion in a SummedAreaTable.
 */
struct SatRegion
{
    const SummedAreaTable* _sat;
    int _x, _y;
    int _w, _h;
    double _sum;

    double getSum() const { return _sum; }

    void create(const int x, const int y, const int w, const int h, const SummedAreaTable* sat)
    {
        _x = x;
        _y = y;
        _w = w;
        _h = h;
        _sat = sat;

        _sum = _sat->sum(x, y, x + (w - 1), y, x + (w - 1), y + (h - 1), x, y + (h - 1));
    }

    // MEDIAN CUT Fastest by far
    // median cut split criteria
    bool splitCriteria(const SatRegion& sat) const
    {
        // median cut
        return sat.getSum() * 2.0 >= getSum();
    }

    void split_w(SatRegion& A) const
    {
        for (int w = 1; w <= _w; ++w)
        {
            A.create(_x, _y, w, _h, _sat);

            // if region left has approximately half the energy of the entire thing stahp
            if (splitCriteria(A))
                break;
        }
    }

    void split_h(SatRegion& A) const
    {
        for (int h = 1; h <= _h; ++h)
        {
            A.create(_x, _y, _w, h, _sat);

            // if region top has approximately half the energy of the entire thing stahp
            if (splitCriteria(A))
                break;
        }
    }

    /**
     * Split region horizontally into subregions A and B.
     */
    void split_w(SatRegion& A, SatRegion& B) const
    {
        split_w(A);
        B.create(_x + (A._w - 1), _y, _w - A._w, _h, _sat);
    }

    /**
     * Split region vertically into subregions A and B.
     */
    void split_h(SatRegion& A, SatRegion& B) const
    {
        split_h(A);
        B.create(_x, _y + (A._h - 1), _w, _h - A._h, _sat);
    }

    double2 centroid() const { return double2(_x + _w * 0.5, _y + _h * 0.5); }
    double areaSize() const { return _w * _h; }
};

struct ExtractLights
{
    int _nbRegions = 0;
    int _nbLights = 0;

    SatRegion* _regions;
    Light* _lights;
    Light* _mainLights;

    Light _mainLight;
    int _nbMainLights = 0;
    int _padding = 0;

    /**
     * Recursively split a region r and append new subregions
     * A and B to regions vector when at an end.
     */
    void splitRecursive(const SatRegion& r, const int n)
    {
        // check: can't split any further?
        if (r._w < 2 || r._h < 2 || n == 0)
        {
            // only now add region
            _regions[_nbRegions++] = r;
            return;
        }

        SatRegion A, B;

        if (r._w > r._h)
            r.split_w(A, B);
        else
            r.split_h(A, B);

        if (A._h > 2 && A._w > 2)
        {
            splitRecursive(A, n - 1);
        }

        if (B._h > 2 && B._w > 2)
        {
            splitRecursive(B, n - 1);
        }
    }

    /**
     * The median cut algorithm Or Variance Minimisation
     *
     * img - Summed area table of an image
     * n - number of subdivision, yields 2^n cuts
     * regions - an empty vector that gets filled with generated regions
     */
    void medianVarianceCut(const SummedAreaTable& img, const int n)
    {
        _nbRegions = 0;

        // insert entire image as start region
        SatRegion region;
        region.create(0, 0, img.width(), img.height(), &img);

        // recursively split into subregions
        splitRecursive(region, n);
    }

    void createLightsFromRegions(const double maxLum, const Image& image, const SummedAreaTable& lumSat)
    {
        const double weigth = lumSat.getWeightAccumulation();
        int width = lumSat.width();
        int height = lumSat.height();
        const int imgSize = width * height;
        double weight = (4.0 * M_PI) / ((double)(imgSize));

        const float3* rgb = static_cast<const float3*>(&image.getPixel(0, 0));

        // convert region into lights
        for (int regionIndex = 0; regionIndex < _nbRegions; regionIndex++)
        {
            const SatRegion& region = _regions[regionIndex];

            Light& light = _lights[_nbLights++];

            // init values
            light._merged = false;
            light._mergedNum = 0;

            light._position[0] = region._x;
            light._position[1] = region._y;
            light._size[0] = region._w;
            light._size[1] = region._h;

            // set light at centroid
            light._centroidPosition = region.centroid();

            // light area Size
            light._areaSize = region.areaSize();

            // compute area values, as SAT introduce precision errors
            // due to high sum values against small data values
            // we use here less error inducing computations

            int x = static_cast<int>(light._position[0]);
            int y = static_cast<int>(light._position[1]);
            double3 colorSum(0.0, 0.0, 0.0);
            double lumSum = 0.0;

            int lightHeight = (int)light._size[1];
            int lightWidth = (int)light._size[0];
            for (int y1 = y; y1 < y + lightHeight; ++y1)
            {
                const double posY = ((double)y1 + 1.0) / (double)(height + 1.0);
                const double solidAngle = cos(M_PI * (posY - 0.5)) * weight;

                double3 color;
                for (int x1 = x; x1 < x + lightWidth; ++x1)
                {
                    const int i = (x1 + (y1 * width));

                    float3ToDouble3(color.ptr(), rgb[i].ptr());
                    lumSum += luminance(color) * solidAngle;
                    colorSum += color;
                }
            }

            // normalize
            lumSum *= (4.0 * M_PI) / weigth;
            light._sum = lumSum;

            // Colors
            light._colorAverage = colorSum * (1.0 / light._areaSize);
            light._lumAverage = lumSum / light._areaSize;

            // make all value 0..1 now
            light._position[0] = static_cast<double>(light._position[0]) / (double)width;
            light._position[1] = static_cast<double>(light._position[1]) / (double)height;
            light._size[0] = static_cast<double>(light._size[0]) / (double)width;
            light._size[1] = static_cast<double>(light._size[1]) / (double)height;
            light._areaSize = light._size[0] * light._size[1];

            light._centroidPosition[0] = light._centroidPosition[0] / (double)width;
            light._centroidPosition[1] = light._centroidPosition[1] / (double)height;

            // if value out of bounds
            light._error = light._sum > maxLum;
        }
    }

    /*
     *  Merge a light into another and store a copy inside the parent
     */
    void _mergeLight(Light& lightParent, Light& lightChild)
    {
        // exclude from next merges
        lightChild._merged = true;

        const double x = lightParent._position[0];
        const double y = lightParent._position[1];
        const double w = lightParent._size[0];
        const double h = lightParent._size[1];

        lightParent._position[0] = fmin(x, lightChild._position[0]);
        lightParent._position[1] = fmin(y, lightChild._position[1]);

        lightParent._size[0] = fmax(x + w, (lightChild._position[0] + lightChild._size[0])) - lightParent._position[0];
        lightParent._size[1] = fmax(y + h, (lightChild._position[1] + lightChild._size[1])) - lightParent._position[1];

        // light is bigger, better candidate to main light
        lightParent._mergedNum++;

        lightParent._sum += lightChild._sum;

        double newAreaSize = lightParent._areaSize + lightChild._areaSize;
        double ratioParent = lightParent._areaSize / newAreaSize;
        double ratioChild = lightChild._areaSize / newAreaSize;

        lightParent._colorAverage = lightParent._colorAverage * ratioParent + lightChild._colorAverage * ratioChild;

        lightParent._areaSize = newAreaSize;
        lightParent._lumAverage = lightParent._sum / newAreaSize;
    }

    /**
     * Merge small area light neighbour with small area light neighbours
     */
    int mergeLights(Light& mainLight, double lengthSizeMax, double degreeMerge)
    {
        // discard or keep Light too near an current light
        const double border = degreeMerge * M_PI / 360.0;
        int numMergedLightTotal = 0;

        // for each light we try to merge with all other intersecting lights
        // that are in the same neighborhood of the sorted list of lights
        // where neighbors are of near same values
        for (int lightIndex = 0; lightIndex < _nbLights; lightIndex++)
        {
            Light light = _lights[lightIndex];

            // already merged in a previous light
            // we do nothing
            if (light._merged)
                continue;

            // double x1 = lightCurrent._x - border;
            // double y1 = lightCurrent._y - border;
            // double x2 = x1 + lightCurrent._w + border;
            // double y2 = y1 + lightCurrent._h + border;
            double2 pos1 = light._position - double2(border, border);
            double2 pos2 = pos1 + light._size + double2(border, border);

            int numMergedLight;

            do
            {
                numMergedLight = 0;

                // could start at current light
                // lights is sorted by areasize from small to big
                for (int i = 0; i < _nbLights; i++)
                {
                    Light& lightToMerge = _lights[i];

                    // ignore already merged into another and itself
                    if (lightToMerge._merged || i == lightIndex)
                        continue;

                    // if merged do new size will be problematic
                    double2 newPosition(fmin(light._position[0], lightToMerge._position[0]),
                                        fmin(light._position[1], lightToMerge._position[1]));

                    double newParentSizeW =
                        fmax(light._position[0] + light._size[0], (lightToMerge._position[0] + lightToMerge._size[0])) -
                        newPosition[0];
                    if (lengthSizeMax < newParentSizeW)
                        continue;

                    double newParentSizeH =
                        fmax(light._position[1] + light._size[1], (lightToMerge._position[1] + lightToMerge._size[1])) -
                        newPosition[1];

                    if (lengthSizeMax < newParentSizeH)
                        continue;

                    bool intersect2D = !(lightToMerge._position[1] > pos2[1] ||
                                         lightToMerge._position[1] + lightToMerge._size[1] < pos1[1] ||
                                         lightToMerge._position[0] > pos2[0] ||
                                         lightToMerge._position[0] + lightToMerge._size[0] < pos1[0]);
                    //  share borders
                    if (intersect2D)
                    {
                        _mergeLight(light, lightToMerge);

                        pos1 = light._position - double2(border, border);
                        // X1 = light._x - border;
                        // y1 = light._y - border;

                        pos2 = pos1 + light._size + double2(border, border);
                        // x2 = x1 + light._w + border;
                        // y2 = y1 + light._h + border;

                        numMergedLight++;
                        numMergedLightTotal++;
                    }
                }

                // if we're merging we're changing borders
                // means we have new neighbours
                // or light now included inside our area
            } while (numMergedLight > 0);

            if (light._mergedNum > 0)
            {
                _mainLights[_nbMainLights++] = light;
            }
        }

        // count merged light
        numMergedLightTotal = 0;
        for (int i = 0; i < _nbMainLights; i++)
        {
            if (!_mainLights[i]._merged)
                numMergedLightTotal += _mainLights[i]._mergedNum;
        }

        // fill new array with light that wasn't merged at all
        for (int i = 0; i < _nbLights; i++)
        {
            Light& lightCurrent = _lights[i];
            // add remaining non merged lights
            if (!lightCurrent._merged && lightCurrent._mergedNum == 0)
            {

                lightCurrent._lumAverage = lightCurrent._sum / lightCurrent._areaSize;
                _mainLights[_nbMainLights++] = lightCurrent;
            }
        }

        // sort By sum now (changed the sort Criteria during merge)
        // biggest Sum first
        qsort(_mainLights, _nbMainLights, sizeof(Light), lightInvCmpSum);

        bool foundLight = false;
        for (int i = 0; i < _nbMainLights; i++)
        {
            if (_mainLights[i]._centroidPosition[1] >= 0.5)
                continue;

            mainLight = _mainLights[i];
            foundLight = true;
            break;
        }

        return foundLight ? 0 : 1;
    }

    int extractMainLight(Light& mainLight, const Image& image, int numCuts = 8)
    {
        const double ratioLuminanceLight = 0.5;

        int nbElements = 2 << numCuts;
        _regions = new SatRegion[nbElements];
        _lights = new Light[nbElements];
        _mainLights = new Light[nbElements];

        ////////////////////////////////////////////////
        // create summed area table of luminance image
        SummedAreaTable sat;
        sat.create(image);

        ////////////////////////////////////////////////
        // apply cut algorithm
        medianVarianceCut(sat, numCuts); // max 2^n cuts

        if (!_nbRegions)
        {
            printf("Cannot cut into light regions\n");
            return 1;
        }

        ////////////////////////////////////////////////
        // create Lights from regions

        // convert absolute input parameters
        // to relative to environment at hand value.
        // From ratio to pixel squared area
        /// Light Max luminance in percentage
        double luminanceSum = sat.getSumLuminance();
        const double luminanceMaxLight = ratioLuminanceLight * luminanceSum;

        // And he saw that light was good, and separated light from darkness
        createLightsFromRegions(luminanceMaxLight, image, sat);

        // sort lights
        // the smaller, the more powerful luminance
        qsort(_lights, _nbLights, sizeof(Light), lightCmpAreaSize);

        const double degreeMerge = 35.0;

        // Light Area Size under which we merge
        // default to size of the median region size
        // if lots of small lights => give small area
        // if lots of big lights => give big area
        // const uint mergeindexPos =  (lights.size() * 25) / 100;
        // const double mergeAreaSize = lights[mergeindexPos]._areaSize;

        double ratioLengthSizeMax = 0.08;
        int errorMergeLight = mergeLights(mainLight, ratioLengthSizeMax, degreeMerge);

        sat.free();
        delete[] _regions;
        delete[] _lights;
        delete[] _mainLights;

        // normalise light color
        mainLight._colorAverage.normalize();

        return errorMergeLight;
    }
};

int extractMainLight(Light& light, const Image& equirect)
{
    ExtractLights extractLight;
    return extractLight.extractMainLight(light, equirect);
}
