#pragma once

#include <stdlib.h>
#include <float.h>

// template

struct double2
{
    double _v[2];

    double2()
        : _v{0, 0}
    {}
    double2(double x, double y)
        : _v{x, y}
    {}

    inline double& operator[](int i) { return _v[i]; }
    inline const double& operator[](int i) const { return _v[i]; }

    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline double2& operator+=(const double2& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }

    /** Multiply by scalar. */
    inline const double2 operator*(float rhs) const { return double2(_v[0] * rhs, _v[1] * rhs); }
};

template <typename T> struct Vec3Type
{
    /** Data type of vector components.*/
    typedef T value_type;

    value_type _v[3];

    Vec3Type()
        : _v{0, 0, 0}
    {}
    Vec3Type(value_type x, value_type y, value_type z)
        : _v{x, y, z}
    {}

    inline value_type& operator[](int i) { return _v[i]; }
    inline const value_type& operator[](int i) const { return _v[i]; }

    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline Vec3Type& operator+=(const Vec3Type& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        _v[2] += rhs._v[2];
        return *this;
    }

    /** Multiply by scalar. */
    inline Vec3Type& operator*=(const value_type& rhs)
    {
        _v[0] *= rhs;
        _v[1] *= rhs;
        _v[2] *= rhs;
        return *this;
    }

    /** Multiply by scalar. */
    inline const Vec3Type operator*(value_type rhs) const { return Vec3Type(_v[0] * rhs, _v[1] * rhs, _v[2] * rhs); }
};

typedef Vec3Type<float> float3;
typedef Vec3Type<double> double3;

template <typename T> struct Vec4Type
{
  public:
    /** Data type of vector components.*/
    typedef T value_type;

    /** Vec member variable. */
    value_type _v[4];

    // Methods are defined here so that they are implicitly inlined

    /** Constructor that sets all components of the vector to zero */
    Vec4Type()
        : _v{0, 0, 0, 0}
    {}

    Vec4Type(value_type x, value_type y, value_type z, value_type w)
        : _v{x, y, z, w}
    {}

    Vec4Type(const Vec3Type<T>& v3, value_type w)
        : _v{v3[0], v3[1], v3[2], w}
    {}

    inline value_type* ptr() { return _v; }
    inline const value_type* ptr() const { return _v; }

    inline value_type& operator[](unsigned int i) { return _v[i]; }
    inline value_type operator[](unsigned int i) const { return _v[i]; }

    /** Returns true if all components have values that are not NaN. */
    inline bool valid() const { return !isNaN(); }
    /** Returns true if at least one component has value NaN. */
    inline bool isNaN() const { return isnan(_v[0]) || isnan(_v[1]) || isnan(_v[2]) || isnan(_v[3]); }

    /** Dot product. */
    inline value_type operator*(const Vec4Type& rhs) const
    {
        return _v[0] * rhs._v[0] + _v[1] * rhs._v[1] + _v[2] * rhs._v[2] + _v[3] * rhs._v[3];
    }

    /** Multiply by scalar. */
    inline Vec4Type operator*(value_type rhs) const
    {
        return Vec4Type(_v[0] * rhs, _v[1] * rhs, _v[2] * rhs, _v[3] * rhs);
    }

    /** Unary multiply by scalar. */
    inline Vec4Type& operator*=(value_type rhs)
    {
        _v[0] *= rhs;
        _v[1] *= rhs;
        _v[2] *= rhs;
        _v[3] *= rhs;
        return *this;
    }

    /** Divide by scalar. */
    inline Vec4Type operator/(value_type rhs) const
    {
        return Vec4Type(_v[0] / rhs, _v[1] / rhs, _v[2] / rhs, _v[3] / rhs);
    }

    /** Unary divide by scalar. */
    inline Vec4Type& operator/=(value_type rhs)
    {
        _v[0] /= rhs;
        _v[1] /= rhs;
        _v[2] /= rhs;
        _v[3] /= rhs;
        return *this;
    }

    /** Binary vector add. */
    inline Vec4Type operator+(const Vec4Type& rhs) const
    {
        return Vec4Type(_v[0] + rhs._v[0], _v[1] + rhs._v[1], _v[2] + rhs._v[2], _v[3] + rhs._v[3]);
    }

    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline Vec4Type& operator+=(const Vec4Type& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        _v[2] += rhs._v[2];
        _v[3] += rhs._v[3];
        return *this;
    }

    /** Binary vector subtract. */
    inline Vec4Type operator-(const Vec4Type& rhs) const
    {
        return Vec4Type(_v[0] - rhs._v[0], _v[1] - rhs._v[1], _v[2] - rhs._v[2], _v[3] - rhs._v[3]);
    }

    /** Unary vector subtract. */
    inline Vec4Type& operator-=(const Vec4Type& rhs)
    {
        _v[0] -= rhs._v[0];
        _v[1] -= rhs._v[1];
        _v[2] -= rhs._v[2];
        _v[3] -= rhs._v[3];
        return *this;
    }

    /** Negation operator. Returns the negative of the Vec4Type. */
    inline const Vec4Type operator-() const { return Vec4Type(-_v[0], -_v[1], -_v[2], -_v[3]); }

    /** Length of the vector = sqrt( vec . vec ) */
    inline value_type length() const { return sqrt(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2] + _v[3] * _v[3]); }

    /** Length squared of the vector = vec . vec */
    inline value_type length2() const { return _v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2] + _v[3] * _v[3]; }

    /** Normalize the vector so that it has length unity.
     * Returns the previous length of the vector.
     */
    inline value_type normalize()
    {
        value_type norm = Vec4Type::length();
        if (norm > 0.0)
        {
            value_type inv = 1.0 / norm;
            _v[0] *= inv;
            _v[1] *= inv;
            _v[2] *= inv;
            _v[3] *= inv;
        }
        return (norm);
    }

}; // end of class Vec4Type

typedef Vec4Type<float> float4;
typedef Vec4Type<double> double4;

inline double2 hammersley(uint32_t i, double inverseSampleCount)
{

    uint32_t bits = i;
    bits = (bits << 16) | (bits >> 16);
    bits = ((bits & 0x55555555) << 1) | ((bits & 0xAAAAAAAA) >> 1);
    bits = ((bits & 0x33333333) << 2) | ((bits & 0xCCCCCCCC) >> 2);
    bits = ((bits & 0x0F0F0F0F) << 4) | ((bits & 0xF0F0F0F0) >> 4);
    bits = ((bits & 0x00FF00FF) << 8) | ((bits & 0xFF00FF00) >> 8);
    double y = double(bits) * 2.3283064365386963e-10; // / 0x100000000

    return double2(double(i) * inverseSampleCount, y);
}
