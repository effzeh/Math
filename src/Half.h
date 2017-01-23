#pragma once

#define LIB_MATH_HALF_H 1
#define LIB_MATH_HALF_USE_F16C_INSTRINSICS 1

#include <cstdint>
#if LIB_MATH_HALF_USE_F16C_INSTRINSICS
#include <immintrin.h>
#endif

namespace math {

//==============================================================================
// Conversions
//==============================================================================

// CUDA: __half2float
// CUDA: __float2half_rn

union F16Bits
{
    uint16_t u;
#ifdef __GNUC__
    __extension__
#endif
    struct // assumes little-endian!
    {
        uint16_t mantissa : 10;
        uint16_t exponent : 5;
        uint16_t sign : 1;
    };
};

union F32Bits
{
    uint32_t u;
#ifdef __GNUC__
    __extension__
#endif
    struct // assumes little-endian!
    {
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    };
    float f;
};

inline float Float16ToFloat32(uint16_t value)
{
#if LIB_MATH_HALF_USE_F16C_INSTRINSICS
    auto xmm0 = _mm_cvtsi32_si128(value);
    auto xmm1 = _mm_cvtph_ps(xmm0);

    return _mm_cvtss_f32(xmm1);
#else
    // From:
    // https://gist.github.com/rygorous/2144712

    static const F32Bits magic = { 126 << 23 };

    F16Bits h = { value };
    F32Bits f;

    if (h.exponent == 0) // Zero / Denormal
    {
        f.u  = magic.u + h.mantissa;
        f.f -= magic.f;
    }
    else
    {
        f.mantissa = h.mantissa << 13;

        if (h.exponent == 0x1f) // Inf/NaN
        {
            f.exponent = 255;
            if (h.mantissa) // NaN -> qNaN
            {
                f.mantissa |= 1 << 22;
            }
        }
        else
        {
            f.exponent = 127 - 15 + h.exponent;
        }
    }

    f.sign = h.sign;

    return f.f;
#endif
}

inline uint16_t Float32ToFloat16(float value)
{
#if LIB_MATH_HALF_USE_F16C_INSTRINSICS
    auto xmm0 = _mm_load_ss(&value);
    auto xmm1 = _mm_cvtps_ph(xmm0, 0);

    return static_cast<uint16_t>(_mm_cvtsi128_si32(xmm1));
#else
    // From:
    // https://gist.github.com/rygorous/2156668

    static const F32Bits inf32 = { 255 << 23 };
    static const F32Bits inf16 = {  31 << 23 };
    static const F32Bits magic = {  15 << 23 };

    F32Bits f;
    F16Bits h;

    f.f = value;

    uint32_t sign = f.u & 0x80000000;

    f.u ^= sign;

    // All the integer compares in this function can be safely compiled
    // into signed compares since all operands are below 0x80000000. Important
    // if you want fast straight SSE2 code (since there's no unsigned PCMPGTD).

    if (f.u >= inf32.u) // Inf or NaN (all exponent bits set)
    {
        h.u = 0x7c00 | ((f.u >> 13) & 0x3FF);
    }
    else // (De)normalized number or zero
    {
        uint32_t round_mask = ~0xfffu;

        f.u &= round_mask;
        f.f *= magic.f;
        f.u -= round_mask;

        if (f.u > inf16.u) // Clamp to signed infinity if overflowed
        {
            f.u = inf16.u;
        }

        h.u = (uint16_t)(f.u >> 13); // Take the bits!
    }

    h.u |= sign >> 16;

    return h.u;
#endif
}

//==============================================================================
// Types
//==============================================================================

struct Half
{
    union
    {
        uint16_t u;
#ifdef __GNUC__
        __extension__
#endif
        struct // FIXME: assumes little-endian!
        {
            uint16_t mantissa : 10;
            uint16_t exponent : 5;
            uint16_t sign : 1;
        };
    };

    Half() = default;

    // Construct from float
    /*implicit*/ Half(float f) : u(Float32ToFloat16(f))
    {
    }

    // Convert to float
    /*implicit*/ operator float() const
    {
        return Float16ToFloat32(u);
    }

    // Construct from bits
    static Half fromBits(uint16_t u)
    {
        Half h;
        h.u = u;
        return h;
    }

    // Construct from bits
    static Half fromBits(uint16_t mantissa, uint16_t exponent, uint16_t sign)
    {
        Half h;

        h.mantissa = mantissa;
        h.exponent = exponent;
        h.sign = sign;

        return h;
    }
};

template <class F> Half& operator +=(Half& lhs, F const& rhs) { return lhs = lhs + rhs; }
template <class F> Half& operator -=(Half& lhs, F const& rhs) { return lhs = lhs - rhs; }
template <class F> Half& operator *=(Half& lhs, F const& rhs) { return lhs = lhs * rhs; }
template <class F> Half& operator /=(Half& lhs, F const& rhs) { return lhs = lhs / rhs; }

} // namespace math
