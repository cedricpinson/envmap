#pragma once

#include <float.h>
#include <math.h>
#include <stdlib.h>

// template

template <typename T> struct Vec2Type
{
    /** Data type of vector components.*/
    typedef T value_type;

    /** Vec member variable. */
    value_type _v[2];

    /** Constructor that sets all components of the vector to zero */
    Vec2Type()
    {
        _v[0] = 0.0;
        _v[1] = 0.0;
    }
    Vec2Type(value_type x, value_type y)
    {
        _v[0] = x;
        _v[1] = y;
    }

    inline bool operator==(const Vec2Type& v) const { return _v[0] == v._v[0] && _v[1] == v._v[1]; }

    inline bool operator!=(const Vec2Type& v) const { return _v[0] != v._v[0] || _v[1] != v._v[1]; }

    inline bool operator<(const Vec2Type& v) const
    {
        if (_v[0] < v._v[0])
            return true;
        else if (_v[0] > v._v[0])
            return false;
        else
            return (_v[1] < v._v[1]);
    }

    inline value_type* ptr() { return _v; }
    inline const value_type* ptr() const { return _v; }

    inline void set(value_type x, value_type y)
    {
        _v[0] = x;
        _v[1] = y;
    }

    inline value_type& operator[](int i) { return _v[i]; }
    inline value_type operator[](int i) const { return _v[i]; }

    inline value_type& x() { return _v[0]; }
    inline value_type& y() { return _v[1]; }

    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }

    /** Dot product. */
    inline value_type operator*(const Vec2Type& rhs) const { return _v[0] * rhs._v[0] + _v[1] * rhs._v[1]; }

    /** Multiply by scalar. */
    inline const Vec2Type operator*(value_type rhs) const { return Vec2Type(_v[0] * rhs, _v[1] * rhs); }

    /** Unary multiply by scalar. */
    inline Vec2Type& operator*=(value_type rhs)
    {
        _v[0] *= rhs;
        _v[1] *= rhs;
        return *this;
    }

    /** Divide by scalar. */
    inline const Vec2Type operator/(value_type rhs) const { return Vec2Type(_v[0] / rhs, _v[1] / rhs); }

    /** Unary divide by scalar. */
    inline Vec2Type& operator/=(value_type rhs)
    {
        _v[0] /= rhs;
        _v[1] /= rhs;
        return *this;
    }

    /** Binary vector add. */
    inline const Vec2Type operator+(const Vec2Type& rhs) const
    {
        return Vec2Type(_v[0] + rhs._v[0], _v[1] + rhs._v[1]);
    }

    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline Vec2Type& operator+=(const Vec2Type& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }

    /** Binary vector subtract. */
    inline const Vec2Type operator-(const Vec2Type& rhs) const
    {
        return Vec2Type(_v[0] - rhs._v[0], _v[1] - rhs._v[1]);
    }

    /** Unary vector subtract. */
    inline Vec2Type& operator-=(const Vec2Type& rhs)
    {
        _v[0] -= rhs._v[0];
        _v[1] -= rhs._v[1];
        return *this;
    }

    /** Negation operator. Returns the negative of the Vec2Type. */
    inline const Vec2Type operator-() const { return Vec2Type(-_v[0], -_v[1]); }

    /** Length of the vector = sqrt( vec . vec ) */
    inline value_type length() const { return sqrtf(_v[0] * _v[0] + _v[1] * _v[1]); }

    /** Length squared of the vector = vec . vec */
    inline value_type length2(void) const { return _v[0] * _v[0] + _v[1] * _v[1]; }

    /** Normalize the vector so that it has length unity.
     * Returns the previous length of the vector.
     */
    inline value_type normalize()
    {
        value_type norm = Vec2Type::length();
        if (norm > 0.0)
        {
            value_type inv = 1.0f / norm;
            _v[0] *= inv;
            _v[1] *= inv;
        }
        return (norm);
    }
};

typedef Vec2Type<float> float2;
typedef Vec2Type<double> double2;

template <typename T> struct Vec3Type
{
    /** Data type of vector components.*/
    typedef T value_type;

    value_type _v[3];

    /** Constructor that sets all components of the vector to zero */
    Vec3Type()
        : _v{0, 0, 0}
    {}

    Vec3Type(value_type x, value_type y, value_type z)
        : _v{x, y, z}
    {}

    inline bool operator==(const Vec3Type& v) const { return _v[0] == v._v[0] && _v[1] == v._v[1] && _v[2] == v._v[2]; }

    inline bool operator!=(const Vec3Type& v) const { return _v[0] != v._v[0] || _v[1] != v._v[1] || _v[2] != v._v[2]; }

    inline bool operator<(const Vec3Type& v) const
    {
        if (_v[0] < v._v[0])
            return true;
        else if (_v[0] > v._v[0])
            return false;
        else if (_v[1] < v._v[1])
            return true;
        else if (_v[1] > v._v[1])
            return false;
        else
            return (_v[2] < v._v[2]);
    }

    inline value_type* ptr() { return _v; }
    inline const value_type* ptr() const { return _v; }

    Vec3Type<double> toDouble() const { return Vec3Type<double>((double)_v[0], (double)_v[1], (double)_v[2]); }

    inline void set(value_type x, value_type y, value_type z)
    {
        _v[0] = x;
        _v[1] = y;
        _v[2] = z;
    }

    inline void set(const Vec3Type& rhs)
    {
        _v[0] = rhs._v[0];
        _v[1] = rhs._v[1];
        _v[2] = rhs._v[2];
    }

    inline Vec3Type max(value_type m) { return Vec3Type(fmax(_v[0], m), fmax(_v[1], m), fmax(_v[2], m)); }

    inline value_type& operator[](int i) { return _v[i]; }
    inline value_type operator[](int i) const { return _v[i]; }

    inline value_type& x() { return _v[0]; }
    inline value_type& y() { return _v[1]; }
    inline value_type& z() { return _v[2]; }

    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }
    inline value_type z() const { return _v[2]; }

    /** Dot product. */
    inline value_type operator*(const Vec3Type& rhs) const
    {
        return _v[0] * rhs._v[0] + _v[1] * rhs._v[1] + _v[2] * rhs._v[2];
    }

    /** Cross product. */
    inline const Vec3Type operator^(const Vec3Type& rhs) const
    {
        return Vec3Type(_v[1] * rhs._v[2] - _v[2] * rhs._v[1], _v[2] * rhs._v[0] - _v[0] * rhs._v[2],
                        _v[0] * rhs._v[1] - _v[1] * rhs._v[0]);
    }

    /** Multiply by scalar. */
    inline const Vec3Type operator*(value_type rhs) const { return Vec3Type(_v[0] * rhs, _v[1] * rhs, _v[2] * rhs); }

    /** Unary multiply by scalar. */
    inline Vec3Type& operator*=(value_type rhs)
    {
        _v[0] *= rhs;
        _v[1] *= rhs;
        _v[2] *= rhs;
        return *this;
    }

    /** Divide by scalar. */
    inline const Vec3Type operator/(value_type rhs) const { return Vec3Type(_v[0] / rhs, _v[1] / rhs, _v[2] / rhs); }

    /** Unary divide by scalar. */
    inline Vec3Type& operator/=(value_type rhs)
    {
        _v[0] /= rhs;
        _v[1] /= rhs;
        _v[2] /= rhs;
        return *this;
    }

    /** Binary vector add. */
    inline const Vec3Type operator+(const Vec3Type& rhs) const
    {
        return Vec3Type(_v[0] + rhs._v[0], _v[1] + rhs._v[1], _v[2] + rhs._v[2]);
    }

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

    /** Binary vector subtract. */
    inline const Vec3Type operator-(const Vec3Type& rhs) const
    {
        return Vec3Type(_v[0] - rhs._v[0], _v[1] - rhs._v[1], _v[2] - rhs._v[2]);
    }

    /** Unary vector subtract. */
    inline Vec3Type& operator-=(const Vec3Type& rhs)
    {
        _v[0] -= rhs._v[0];
        _v[1] -= rhs._v[1];
        _v[2] -= rhs._v[2];
        return *this;
    }

    /** Negation operator. Returns the negative of the Vec3Type. */
    inline const Vec3Type operator-() const { return Vec3Type(-_v[0], -_v[1], -_v[2]); }

    /** Length of the vector = sqrt( vec . vec ) */
    inline value_type length() const { return sqrtf(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]); }

    /** Length squared of the vector = vec . vec */
    inline value_type length2() const { return _v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]; }

    /** Normalize the vector so that it has length unity.
     * Returns the previous length of the vector.
     */
    inline value_type normalize()
    {
        value_type norm = Vec3Type::length();
        if (norm > 0.0)
        {
            value_type inv = 1.0f / norm;
            _v[0] *= inv;
            _v[1] *= inv;
            _v[2] *= inv;
        }
        return (norm);
    }
}; // end of class Vec3Type

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

    return double2(i * inverseSampleCount, y);
}

inline double D_GGX(double NoH, double alpha)
{
    // use GGX / Trowbridge-Reitz, same as Disney and Unreal 4
    // cf http://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_notes_v2.pdf p3
    // const double PI_INV = 1.0 / M_PI;
    // double tmp = alpha / (NdotH * NdotH * (alpha * alpha - 1.0) + 1.0);
    // return tmp * tmp * PI_INV;
    // from filament note
    // NOTE: (aa-1) == (a-1)(a+1) produces better fp accuracy
    double f = (alpha - 1) * ((alpha + 1) * (NoH * NoH)) + 1.0;
    return (alpha * alpha) / (M_PI * f * f);
}

template <typename T> inline typename T::value_type dot(const T& a, const T& b) { return a * b; }
template <typename T> inline T cross(const T& a, const T& b) { return a ^ b; }
template <typename T> inline T normalize(const T& a) { return a * (1.0f / a.length()); }
template <typename T> inline T lerp(const T& a, const T& b, float t) { return a + ((b - a) * t); }
template <typename T> inline T frac(T v) { return v - floor(v); }
template <typename T> inline T clamp(T v, T minimum, T maximum)
{
    return v < minimum ? minimum : v > maximum ? maximum : v;
}

template <typename T> inline void cross(T* result, const T* a, const T* b)
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}
template <typename T> inline void normalize(T* result, const T* a)
{
    T l = 1.0f / sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    result[0] = a[0] * l;
    result[1] = a[1] * l;
    result[2] = a[2] * l;
}
inline void float3ToDouble3(double* result, const float* a)
{
    result[0] = (double)a[0];
    result[1] = (double)a[1];
    result[2] = (double)a[2];
}
