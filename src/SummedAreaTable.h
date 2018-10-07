#pragma once

struct Image;

/**
 * Luminance Summed Area Table
 * https://en.wikipedia.org/wiki/Summed_area_table
 * Aka Integral Image (for order > 1)
 * Accelerate Sum computation
 * Create a luminance summed area table from an image.
 */
struct SummedAreaTable
{
    int _width, _height;
    double* _sat;

    // error reducing values
    double _weightAccum;
    double _sumLuminance;

    void create(const Image& image);
    void free();

    double getSumLuminance() const { return _sumLuminance; }
    double getWeightAccumulation() const { return _weightAccum; }

    double I(const int x, const int y) const
    {
        if (x < 0 || y < 0)
            return 0.0;
        if (x >= _width || y >= _height)
            return 0.0;
        const int i = y * _width + x;
        return _sat[i];
    }

    int width() const { return _width; }
    int height() const { return _height; }

    /**
     * Returns the sum of a region defined by A,B,C,D.
     *
     * A----B
     * |    |  sum = C+A-B-D
     * D----C
     */
    double sum(const int ax, const int ay, const int bx, const int by, const int cx, const int cy, const int dx,
               const int dy) const
    {
        return I(cx, cy) + I(ax, ay) - I(bx, by) - I(dx, dy);
    }
};
