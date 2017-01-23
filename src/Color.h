#pragma once

#define LIB_MATH_COLOR_CONVERSIONS_H 1

#include "Vector.h"

//
// NOTE:
// If not otherwise stated, RGB here means sRGB (D65) color space.
//

namespace math {

//------------------------------------------------------------------------------
// Linear RGB / Non-linear RGB
//------------------------------------------------------------------------------

// Convert linear RGB to non-linear RGB
inline float LinearRGBToNonlinearRGB(float c)
{
    return (c <= 0.0031308f)
                ? 12.92f * c
                : 1.055f * Pow(c, 1.0f/2.4f) - 0.055f;
}

// Convert linear RGB to non-linear RGB
inline Vec3 LinearRGBToNonlinearRGB(float R, float G, float B) {
    return { LinearRGBToNonlinearRGB(R), LinearRGBToNonlinearRGB(G), LinearRGBToNonlinearRGB(B) };
}

// Convert linear RGB to non-linear RGB
inline Vec3 LinearRGBToNonlinearRGB(Vec3 const& c) {
    return LinearRGBToNonlinearRGB(c.x, c.y, c.z);
}

// Convert non-linear to linear RGB
inline float NonlinearRGBToLinearRGB(float c)
{
    return (c <= 0.04045f)
                ? c / 12.92f
                : Pow((c + 0.055f) / 1.055f, 2.4f);
}

// Convert non-linear to linear RGB
inline Vec3 NonlinearRGBToLinearRGB(float R, float G, float B) {
    return { NonlinearRGBToLinearRGB(R), NonlinearRGBToLinearRGB(G), NonlinearRGBToLinearRGB(B) }; }

// Convert non-linear to linear RGB
inline Vec3 NonlinearRGBToLinearRGB(Vec3 const& c) {
    return NonlinearRGBToLinearRGB(c.x, c.y, c.z);
}

//------------------------------------------------------------------------------
// Linear RGB / HSV
//------------------------------------------------------------------------------

// Convert hue to linear RGB
// H must be in the range [0,1]
inline Vec3 HueToLinearRGB(float h)
{
//  const float t = 6.0f * Saturate(h);
    const float t = 6.0f * h;

    const float R = Saturate( Abs(t - 3.0f) - 1.0f );
    const float G = Saturate( 2.0f - Abs(t - 2.0f) );
    const float B = Saturate( 2.0f - Abs(t - 4.0f) );

    return { R, G, B };
}

// Convert HSV to linear RGB
inline Vec3 HSVToLinearRGB(float H, float S, float V) {
    return ((HueToLinearRGB(H) - 1.0f) * S + 1.0f) * V;
}

// Convert HSV to linear RGB
inline Vec3 HSVToLinearRGB(Vec3 const& c) { return HSVToLinearRGB(c.x, c.y, c.z); }

// Convert linear RGB to HSV
inline Vec3 LinearRGBToHSV(float R, float G, float B)
{
    const float maxRGB = Max(R, G, B);
    const float range = maxRGB - Min(R, G, B);

    float H;
    if (range == 0.0f)
        H = 0.0f;
    else if (maxRGB == R)
        H = (G - B) / range;
    else if (maxRGB == G)
        H = (B - R) / range + 2.0f;
    else if (maxRGB == B)
        H = (R - G) / range + 4.0f;
    else
        H = 0.0f;

    const float S = (maxRGB == 0.0f) ? 0.0f : range / maxRGB;

    const float V = maxRGB;

    return { H, S, V };
}

// Convert linear RGB to HSV
inline Vec3 LinearRGBToHSV(Vec3 const& c) { return LinearRGBToHSV(c.x, c.y, c.z); }

//------------------------------------------------------------------------------
// Linear RGB / XYZ
//------------------------------------------------------------------------------

// Convert linear RGB to XYZ
inline Vec3 LinearRGBToXYZ(float R, float G, float B)
{
    const float X = 0.41238656f * R + 0.35759149f * G + 0.18045049f * B;
    const float Y = 0.21263682f * R + 0.71518298f * G + 0.07218020f * B;
    const float Z = 0.01933062f * R + 0.11919716f * G + 0.95037259f * B;

    return { X, Y, Z };
}

// Convert linear RGB to XYZ
inline Vec3 LinearRGBToXYZ(Vec3 const& c) { return LinearRGBToXYZ(c.x, c.y, c.z); }

// Returns the luminance (i.e. the Y in XYZ)
inline float LinearRGBToLuminance(float R, float G, float B) {
    return 0.21263682f * R + 0.71518298f * G + 0.07218020f * B;
}

// Returns the luminance (i.e. the Y in XYZ)
inline float LinearRGBToLuminance(Vec3 const& c) { return LinearRGBToLuminance(c.x, c.y, c.z); }

// Convert XYZ to linear RGB
inline Vec3 XYZToLinearRGB(float X, float Y, float Z)
{
    const float R =  3.24100323f * X - 1.53739897f * Y - 0.49861588f * Z;
    const float G = -0.96922425f * X + 1.87592998f * Y + 0.04155423f * Z;
    const float B =  0.05563942f * X - 0.20401121f * Y + 1.05714898f * Z;

    return { R, G, B };
}

// Convert XYZ to linear RGB
inline Vec3 XYZToLinearRGB(Vec3 const& c) { return XYZToLinearRGB(c.x, c.y, c.z); }

//------------------------------------------------------------------------------
// Spectrum to XYZ
//------------------------------------------------------------------------------

template <class Emissive, class Observer>
Vec3 SpectrumToXYZ(Emissive const& spd, Observer const& obs)
{
    Vec3 XYZ(0.0f);

    for (int l = 360; l <= 830; ++l)
        XYZ += spd(l) * obs(l); // * dl...

#if 0
    if (XYZ.y == 0.0f)
        return Vec3(0.0f);

    return { XYZ.x / XYZ.y, 1.0f, XYZ.z / XYZ.y };
#else
    return XYZ / XYZ.y;
#endif
}

template <class Reflective, class Illuminant, class Observer>
Vec3 SpectrumToXYZ(Reflective const& spd, Illuminant const& ill, Observer const& obs)
{
    Vec3 XYZ(0.0f);
    float norm = 0.0f;

    for (int l = 360; l <= 830; ++l)
    {
        auto d = ill(l) * obs(l); // * dl...
        XYZ += spd(l) * d;
        norm += d.y;
    }

#if 0
    if (norm == 0.0f)
        return Vec3(0.0f);
#endif

    return XYZ / norm;
}

//------------------------------------------------------------------------------
// Correlated color temperature (CCT)
//------------------------------------------------------------------------------

// Convert correlated color temperature (CCT) to xy.
// Temperature T is in Kelvin.
inline Vec2 TemperatureToXY(float t)
{
#if 0
    // Exact.

    const Vec3 XYZ = SpectrumToXYZ(Blackbody(t), CIEObserver());

    const float x = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
    const float y = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);

    return {x, y};

//  return XYZ.xy() / ReduceAdd(XYZ);
#endif
#if 0
    // Used in the CIE D-series computation

    t = Clamp(t, 4000.0f, 25000.0f);

    const float x = (t <= 7000.0f)
        ? (0.244063f + 0.09911e+3f / t + 2.9678e+6f / (t*t) - 4.6070e+9f / (t*t*t))
        : (0.237040f + 0.24748e+3f / t + 1.9018e+6f / (t*t) - 2.0064e+9f / (t*t*t));

    const float y = -3.0f * (x*x) + 2.87f * x - 0.275f;

    return {x, y};
#endif
#if 1
    // https://en.wikipedia.org/wiki/Planckian_locus#Approximation

    t = Clamp(t, 1000.0f, 15000.0f);

    // Approximate Planckian locus in CIE 1960 UCS
    const float u = (0.860117757f + 1.54118254e-4f * t + 1.28641212e-7f * (t*t))
                  / (1.000000000f + 8.42420235e-4f * t + 7.08145163e-7f * (t*t));
    const float v = (0.317398726f + 4.22806245e-5f * t + 4.20481691e-8f * (t*t))
                  / (1.000000000f - 2.89741816e-5f * t + 1.61456053e-7f * (t*t));

    // Convert to xyY
    const float x = 3.0f * u / (2.0f * u - 8.0f * v + 4.0f);
    const float y = 2.0f * v / (2.0f * u - 8.0f * v + 4.0f);

    return {x, y};
#endif
#if 0
    // Kim et al. 2002
    // https://en.wikipedia.org/wiki/Planckian_locus#Approximation

    t = Clamp(t, 1667.0f, 25000.0f);

    float x;
    float y;

    if (t < 4000.0f)
        x = -0.2661239e+9f / (t*t*t) - 0.2343580e+6f / (t*t) + 0.8776956e+3f / t + 0.179910f;
    else
        x = -3.0258469e+9f / (t*t*t) + 2.1070379e+6f / (t*t) + 0.2226347e+3f / t + 0.240390f;

    if (t < 2222.0f)
        y = -1.1063814f * (x*x*x) - 1.34811020f * (x*x) + 2.18555832f * x - 0.20219683f;
    else if (t < 4000.0f)
        y = -0.9549476f * (x*x*x) - 1.37418593f * (x*x) + 2.09137015f * x - 0.16748867f;
    else
        y =  3.0817580f * (x*x*x) - 5.87338670f * (x*x) + 3.75112997f * x - 0.37001483f;

    return {x, y};
#endif
}

// Convert correlated color temperature (CCT) to linear RGB.
// Temperature T is in Kelvin and should be in the range [1000K, 15000K].
inline Vec3 TemperatureToLinearRGB(float t, float Y = 1.0f)
{
    const Vec2 xy = TemperatureToXY(t);

    Vec3 XYZ;
    XYZ.x = Y * xy.x / xy.y;
    XYZ.y = Y;
    XYZ.z = Y * (1.0f - xy.x - xy.y) / xy.y;

    return XYZToLinearRGB(XYZ);
}

inline Vec3 HeatmapColor(float t)
{
    const float H = Saturate(4.0f/6.0f * (1.0f - t));
    const float s = 0.5f + 0.5f * t;

    return s * HueToLinearRGB(H);
}

//------------------------------------------------------------------------------
// XYZ / CIE Lab/Luv
//------------------------------------------------------------------------------

inline float _cie_f(float x)
{
    const float Threshold = 0.00885645167903563081717167575546f; // (6/29)^3
    const float Scale     = 7.787037037037037037037037037037f;   // 1/3 (29/6)^2

    return (x > Threshold)
                ? Pow(x, 1.0f/3.0f)
                : Scale * x + 4.0f/29.0f;
}

inline float _cie_finv(float x)
{
    const float Threshold = 0.20689655172413793103448275862069f; // 6/29
    const float Scale     = 0.12841854934601664684898929845422f; // 3 * (6/29)^2

    return (x > Threshold)
                ? Pow(x, 3.0f)
                : Scale * (x - 4.0f/29.0f);
}

// Convert XYZ to L*a*b*
inline Vec3 XYZToLab(float X, float Y, float Z, float Xn, float Yn, float Zn)
{
    const float fX = _cie_f(X / Xn);
    const float fY = _cie_f(Y / Yn);
    const float fZ = _cie_f(Z / Zn);

    const float L = 116.0f * fY - 16.0f;
    const float a = 500.0f * (fX - fY);
    const float b = 200.0f * (fY - fZ);

    return { L, a, b };
}

// Convert XYZ to L*a*b*
inline Vec3 XYZToLab(Vec3 const& c, Vec3 const& cn) { return XYZToLab(c.x, c.y, c.z, cn.x, cn.y, cn.z); }

// Convert XYZ to L*a*b*
// Using reference white D65 and 2 deg observer
inline Vec3 XYZToLab(float X, float Y, float Z)
{
    const float Xn =  95.046837f;
    const float Yn = 100.000000f;
    const float Zn = 108.885246f;

    return XYZToLab(X, Y, Z, Xn, Yn, Zn);
}

// Convert XYZ to L*a*b*
// Using reference white D65 and 2 deg observer
inline Vec3 XYZToLab(Vec3 const& c) { return XYZToLab(c.x, c.y, c.z); }

// Convert L*a*b* to XYZ
inline Vec3 LabToXYZ(float L, float a, float b, float Xn, float Yn, float Zn)
{
    const float fY = 1.0f/116.0f * (L + 16.0f);
    const float fX = fY + 1.0f/500.0f * a;
    const float fZ = fY - 1.0f/200.0f * b;

    const float X = Xn * _cie_finv(fX);
    const float Y = Yn * _cie_finv(fY);
    const float Z = Zn * _cie_finv(fZ);

    return { X, Y, Z };
}

// Convert L*a*b* to XYZ
inline Vec3 LabToXYZ(Vec3 const& c, Vec3 const& cn) { return LabToXYZ(c.x, c.y, c.z, cn.x, cn.y, cn.z); }

// Convert L*a*b* to XYZ
// Using reference white D65 and 2 deg observer
inline Vec3 LabToXYZ(float L, float a, float b)
{
    const float Xn =  95.046837f;
    const float Yn = 100.000000f;
    const float Zn = 108.885246f;

    return LabToXYZ(L, a, b, Xn, Yn, Zn);
}

// Convert L*a*b* to XYZ
// Using reference white D65 and 2 deg observer
inline Vec3 LabToXYZ(Vec3 const& c) { return LabToXYZ(c.x, c.y, c.z); }

//------------------------------------------------------------------------------
// Color
//------------------------------------------------------------------------------

// RGBA Color. Linear sRGB.
struct Color
{
    float r;
    float g;
    float b;
    float a;

    Color() = default;
    Color(float r, float g, float b, float a = 1.0f) : r(r), g(g), b(b), a(a) {}

    explicit Color(float v, float a = 1.0f)
        : Color(v, v, v, a) {}

    explicit Color(Vec3 const& v, float a = 1.0f)
        : Color(v.x, v.y, v.z, a) {}

    explicit Color(Vec4 const& v)
        : Color(v.x, v.y, v.z, v.w) {}

    float& operator [](int index) { return (&r)[index]; }
    float const& operator [](int index) const { return (&r)[index]; }

    explicit operator Vec3() const { return { r, g, b }; }
    explicit operator Vec4() const { return { r, g, b, a }; }

    // Returns the luminance of this color.
    // NOTE: Assumes linear RGB!
    float Luminance() const { return LinearRGBToLuminance(r, g, b); }

    // Adjust the luminance of the linear RGB value.
    void SetLuminance(float L)
    {
        const float lum = Luminance();
        if (lum > 0.0f)
        {
            const float scale = L / lum;
            r *= scale;
            g *= scale;
            b *= scale;
        }
    }

    // Adjust the luminance of the linear RGB value.
    Color AdjustLuminance(float L) const
    {
        Color clr = *this;
        clr.SetLuminance(L);
        return clr;
    }

    // Assuming the current values are in non-linear sRGB space, returns the linear RGB values.
    Color ToLinearRGB() const { return Color(NonlinearRGBToLinearRGB(r, g, b), a); }

    // Assuming the current values are in linear RGB space, returns the sRGB values.
    Color ToSRGB() const { return Color(LinearRGBToNonlinearRGB(r, g, b), a); }

    // Assuming the current values are in linear RGB space, returns HSV
    Color ToHSV() const { return Color(LinearRGBToHSV(r, g, b), a); }

    // Assuming the current values are in linear RGB space, returns XYZ
    Color ToXYZ() const { return Color(LinearRGBToXYZ(r, g, b), a); }

    static Color FromLinearRGB(float R, float G, float B, float a = 1.0f) { return Color(R, G, B, a); }
    static Color FromSRGB     (float R, float G, float B, float a = 1.0f) { return Color(NonlinearRGBToLinearRGB(R, G, B), a); }
    static Color FromHSV      (float H, float S, float V, float a = 1.0f) { return Color(HSVToLinearRGB(H, S, V), a); }
    static Color FromXYZ      (float X, float Y, float Z, float a = 1.0f) { return Color(XYZToLinearRGB(X, Y, Z), a); }

    static Color FromLinearRGB(Vec3 const& clr, float a = 1.0f) { return Color(clr, a); }
    static Color FromSRGB     (Vec3 const& clr, float a = 1.0f) { return Color(NonlinearRGBToLinearRGB(clr), a); }
    static Color FromHSV      (Vec3 const& clr, float a = 1.0f) { return Color(HSVToLinearRGB(clr), a); }
    static Color FromXYZ      (Vec3 const& clr, float a = 1.0f) { return Color(XYZToLinearRGB(clr), a); }

    // Convert correlated color temperature to linear RGB.
    // Temperature t is in Kelvin and should be in the range [1000K, 15000K].
    static Color FromTemperature(float t, float a = 1.0f) { return Color(TemperatureToLinearRGB(t), a); }
};

//------------------------------------------------------------------------------
// Compare
//------------------------------------------------------------------------------

inline bool operator ==(Color const& c, Color const& d) {
    return c.r == d.r && c.g == d.g && c.b == d.b && c.a == d.a;
}

inline bool operator !=(Color const& c, Color const& d) {
    return !(c == d);
}

//------------------------------------------------------------------------------
// Arithmetic
//------------------------------------------------------------------------------

inline Color operator +(Color const& c, Color const& d) { return { c.r + d.r, c.g + d.g, c.b + d.b, c.a + d.a }; }
inline Color operator -(Color const& c, Color const& d) { return { c.r - d.r, c.g - d.g, c.b - d.b, c.a - d.a }; }
inline Color operator *(Color const& c, Color const& d) { return { c.r * d.r, c.g * d.g, c.b * d.b, c.a * d.a }; }
inline Color operator /(Color const& c, Color const& d) { return { c.r / d.r, c.g / d.g, c.b / d.b, c.a / d.a }; }

inline Color operator +(Color const& c, float d) { return { c.r + d, c.g + d, c.b + d, c.a + d }; }
inline Color operator -(Color const& c, float d) { return { c.r - d, c.g - d, c.b - d, c.a - d }; }
inline Color operator *(Color const& c, float d) { return { c.r * d, c.g * d, c.b * d, c.a * d }; }
inline Color operator /(Color const& c, float d) { return { c.r / d, c.g / d, c.b / d, c.a / d }; }

inline Color operator +(float c, Color const& d) { return { c + d.r, c + d.g, c + d.b, c + d.a }; }
inline Color operator -(float c, Color const& d) { return { c - d.r, c - d.g, c - d.b, c - d.a }; }
inline Color operator *(float c, Color const& d) { return { c * d.r, c * d.g, c * d.b, c * d.a }; }
inline Color operator /(float c, Color const& d) { return { c / d.r, c / d.g, c / d.b, c / d.a }; }

inline Color& operator +=(Color& c, Color const& d) { c = c + d; return c; }
inline Color& operator -=(Color& c, Color const& d) { c = c - d; return c; }
inline Color& operator *=(Color& c, Color const& d) { c = c * d; return c; }
inline Color& operator /=(Color& c, Color const& d) { c = c / d; return c; }

inline Color& operator +=(Color& c, float d) { c = c + d; return c; }
inline Color& operator -=(Color& c, float d) { c = c - d; return c; }
inline Color& operator *=(Color& c, float d) { c = c * d; return c; }
inline Color& operator /=(Color& c, float d) { c = c / d; return c; }

//------------------------------------------------------------------------------
// Misc.
//------------------------------------------------------------------------------

// Clamp color and alpha values to [0,1]
inline Color Saturate(Color const& c)
{
    const float R = Saturate(c.r);
    const float G = Saturate(c.g);
    const float B = Saturate(c.b);
    const float A = Saturate(c.a);

    return { R, G, B, A };
}

// Returns the luminance of this color.
// NOTE: Assumes linear RGB values!
inline float Luminance(Color const& c) { return c.Luminance(); }

// Premultiply: multiply color values by alpha value
inline Color Premultiply(Color const& c)
{
    const float R = c.r * c.a;
    const float G = c.g * c.a;
    const float B = c.b * c.a;

    return { R, G, B, c.a };
}

// Un-premultiply: divide color values by alpha value
inline Color Unpremultiply(Color const& c)
{
    if (c.a <= 0.0f)
        return {0, 0, 0, 0};

    const float R = c.r / c.a;
    const float G = c.g / c.a;
    const float B = c.b / c.a;

    return { R, G, B, c.a };
}

// Alpha blending.
// Performs the OVER-operation: FG OVER BG.
inline Color BlendPremultiplied(Color const& fg, Color const& bg)
{
#if 1
    return fg + (1.0f - fg.a) * bg;
#else
    const float s = 1.0f - fg.a;

    const float R = fg.r + s * bg.r;
    const float G = fg.g + s * bg.g;
    const float B = fg.b + s * bg.b;
    const float A = fg.a + s * bg.a;

    return { R, G, B, A };
#endif
}

// Alpha blending.
// Performs the OVER-operation: FG OVER BG.
inline Color Blend(Color const& fg, Color const& bg)
{
    const float s = 1.0f - fg.a;

    const float A = fg.a + s * bg.a;
    const float R = (fg.a * fg.r + s * bg.a * bg.r) / A;
    const float G = (fg.a * fg.g + s * bg.a * bg.g) / A;
    const float B = (fg.a * fg.b + s * bg.a * bg.b) / A;

    return { R, G, B, A };
}

} // namespace math
