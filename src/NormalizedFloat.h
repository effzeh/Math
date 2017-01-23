#pragma once

#define LIB_MATH_NORMALIZED_FLOAT_H 1

#include "TMath.h"

namespace math {

template <int Bits, class T> struct Unorm;
template <int Bits, class T> struct Snorm;

using Unorm8    = Unorm<8, uint8_t>;
using Unorm16   = Unorm<16, uint16_t>;
using Snorm8    = Snorm<8, int8_t>;
using Snorm16   = Snorm<16, int16_t>;

//------------------------------------------------------------------------------
// Conversions between (unsigned) normalized and floatint-point.
//
// See:
// OpenGL 4.4 (2.3.4.1/2)
//------------------------------------------------------------------------------

template <int Bits>
inline uint32_t FloatToUnorm(float f)
{
    //static_assert(Bits >= 1 && Bits <= 23, "invalid parameter");

    const float s = (1U << Bits) - 1U;
    return static_cast<uint32_t>( Saturate(f) * s + 0.5f );
}

template <int Bits>
inline float UnormToFloat(uint32_t n)
{
    //static_assert(Bits >= 1 && Bits <= 23, "invalid parameter");

    const float s = 1.0f / ((1U << Bits) - 1U);
    return static_cast<float>( n * s );
}

template <int Bits>
inline int32_t FloatToSnorm(float f)
{
    //static_assert(Bits >= 1 && Bits <= 24, "invalid parameter");

    const float s = (1U << (Bits - 1)) - 1U;
    return static_cast<int32_t>( Round(Clamp(f, -1.0f, 1.0f) * s) );
}

template <int Bits>
inline float SnormToFloat(int32_t n)
{
    //static_assert(Bits >= 1 && Bits <= 24, "invalid parameter");

    const float s = 1.0f / ((1U << (Bits - 1)) - 1U);
    return Max(-1.0f, n * s);
}

//------------------------------------------------------------------------------
// Unorm
//
// From:
// https://msdn.microsoft.com/en-us/library/windows/desktop/dd607323(v=vs.85).aspx
//
// Unsigned normalized integer, meaning that for an n-bit number, all 0's means
// 0.0f, and all 1's means 1.0f. A sequence of evenly spaced floating point
// values from 0.0f to 1.0f are represented. E.g. a 2-bit UNORM represents 0.0f,
// 1/3, 2/3, and 1.0f.
//------------------------------------------------------------------------------

template <int Bits, class T>
struct Unorm
{
    T n;

    Unorm() = default;

    // Construct from unsigned normalized int
    explicit Unorm(uint32_t n)
        : n(n)
    {
    }

    // Construct from float
    explicit Unorm(float f)
        : n(FloatToUnorm<Bits>(f))
    {
    }

    // Convert to unsigned normalized int
    explicit operator uint32_t() const
    {
        return n;
    }

    // Convert to float
    explicit operator float() const
    {
        return UnormToFloat<Bits>(n);
    }

    friend bool operator ==(Unorm lhs, Unorm rhs) { return lhs.n == rhs.n; }
    friend bool operator !=(Unorm lhs, Unorm rhs) { return lhs.n != rhs.n; }
    friend bool operator < (Unorm lhs, Unorm rhs) { return lhs.n <  rhs.n; }
    friend bool operator <=(Unorm lhs, Unorm rhs) { return lhs.n <= rhs.n; }
    friend bool operator > (Unorm lhs, Unorm rhs) { return lhs.n >  rhs.n; }
    friend bool operator >=(Unorm lhs, Unorm rhs) { return lhs.n >= rhs.n; }
};

//------------------------------------------------------------------------------
// Snorm
//
// From:
// https://msdn.microsoft.com/en-us/library/windows/desktop/dd607323(v=vs.85).aspx
//
// Signed normalized integer, meaning that for an n-bit 2's complement number,
// the maximum value means 1.0f (e.g. the 5-bit value 01111 maps to 1.0f), and
// the minimum value means -1.0f (e.g. the 5-bit value 10000 maps to -1.0f). In
// addition, the second-minimum number maps to -1.0f (e.g. the 5-bit value 10001
// maps to -1.0f). There are thus two integer representations for -1.0f. There
// is a single representation for 0.0f, and a single representation for 1.0f.
// This results in a set of integer representations for evenly spaced floating
// point values in the range (-1.0f...0.0f), and also a complementary set of
// representations for numbers in the range (0.0f...1.0f).
//------------------------------------------------------------------------------

template <int Bits, class T>
struct Snorm
{
    T n;

    Snorm() = default;

    // Construct from unsigned normalized int
    explicit Snorm(int32_t n)
        : n(n)
    {
    }

    // Construct from float
    explicit Snorm(float f)
        : n(FloatToSnorm<Bits>(f))
    {
    }

    // Convert to unsigned normalized int
    explicit operator int32_t() const
    {
        return n;
    }

    // Convert to float
    explicit operator float() const
    {
        return SnormToFloat<Bits>(n);
    }

    friend bool operator ==(Snorm lhs, Snorm rhs) { return lhs.n == rhs.n; }
    friend bool operator !=(Snorm lhs, Snorm rhs) { return lhs.n != rhs.n; }
    friend bool operator < (Snorm lhs, Snorm rhs) { return lhs.n <  rhs.n; }
    friend bool operator <=(Snorm lhs, Snorm rhs) { return lhs.n <= rhs.n; }
    friend bool operator > (Snorm lhs, Snorm rhs) { return lhs.n >  rhs.n; }
    friend bool operator >=(Snorm lhs, Snorm rhs) { return lhs.n >= rhs.n; }
};

} // namespace math
