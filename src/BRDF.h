#pragma once

#define LIB_MATH_BRDF_H 1

#include "Vector.h"

namespace math {
namespace brdf {

//==============================================================================
// Diffuse
//==============================================================================

inline Vec3 Diffuse_Lambert(Vec3 const& diff)
{
    return diff * float{kInvPi};
}

inline Vec3 Diffuse_Burley(Vec3 const& diff, float roughness, float NoV, float NoL, float VoH)
{
    // Burley,
    //  "Physically-Based Shading at Disney"

    const float FD90 = 0.5f + 2.0f * VoH * VoH * roughness;
    const float f1 = 1.0f + (FD90 - 1.0f) * Pow5(1 - NoV);
    const float f2 = 1.0f + (FD90 - 1.0f) * Pow5(1 - NoL);

    return diff * (float{kInvPi} * f1 * f2);
}

inline Vec3 Diffuse_OrenNayar(Vec3 const& diff, float roughness, float NoV, float NoL, float VoH)
{
    // Gotanda,
    //  "Beyond a Simple Physically Based Blinn-Phong Model in Real-Time"
    // Fujii,
    //  "A tiny improvement of Oren-Nayar reflectance model"

    const float a = roughness * roughness;

    const float rho = 1.0f + 0.5f * roughness;
    //const float rho = Mix(0.5f, 0.97f, roughness);
    const float s = a;
    const float s2 = s * s;

    const float VoL = 2.0f * VoH * VoH - 1;

    const float u = VoL - NoV * NoL;
    const float v = SelectPositive(u, Max(NoL, NoV), 1);

//  const float A = 1 - 0.5f * s2 / (s2 + 0.33f) + 0.17f * s2 / (s2 + 0.13f) * rho;
    const float A = 1 - 0.5f * s2 / (s2 + 0.33f);
    const float B = 0.45f * s2 / (s2 + 0.09f);

    return diff * (float{kInvPi} * (A + B * u / v) * rho);
}

//==============================================================================
// Specular BRDFs from:
// http://graphicrants.blogspot.de/2013/08/specular-brdf-reference.html
//
// Microfacet BRDF
//      = D G F / (4 N.V N.L)
//      = D G' F,
//
// where G' = G / (4 N.V N.L)
//==============================================================================

//------------------------------------------------------------------------------
// Specular D - Normal distribution function
//------------------------------------------------------------------------------

inline float D_BlinnPhong(float roughness, float NoH)
{
    // Blinn,
    //  "Models of light reflection for computer synthesized pictures"

    const float a = roughness * roughness;
    const float a2 = a * a;

    const float n = 2.0f / a2 - 2.0f;
    const float p = Pow(NoH, n);

    return p * float{kInvPi} / a2; // p * (n + 2.0f) * kInvTwoPi;
}

inline float D_Beckmann(float roughness, float NoH)
{
    // Beckmann,
    //  "The scattering of electromagnetic waves from rough surfaces"

    const float a = roughness * roughness;
    const float a2 = a * a;

    const float d2 = NoH * NoH;

    return Exp((d2 - 1.0f) / (a2 * d2)) / (float{kPi} * a2 * d2 * d2);
}

inline float D_GGX(float roughness, float NoH)
{
    // Walter et al.,
    //  "Microfacet models for refraction through rough surfaces"
    //
    // GGX (Trowbridge-Reitz)

    const float a = roughness * roughness;
    const float a2 = a * a;

    const float d = NoH * (NoH * a2 - NoH) + 1;

    return a2 / (float{kPi} * d * d);
}

//------------------------------------------------------------------------------
// Specular G - Geometric shadowing
//
// Visibility term Vis = G / (4 N.V N.L)
//------------------------------------------------------------------------------

inline float Vis_Implicit()
{
    return 0.25;
}

inline float Vis_Neumann(float NoV, float NoL)
{
    // Neumann et al.,
    //  "Compact metallic reflectance models"

    return 0.25f / Max(NoL, NoV);
}

inline float Vis_Kelemen(float VoH)
{
    // Kelemen,
    //  "A microfacet based coupled specular-matte brdf model with importance sampling"

    return 0.25f / (VoH * VoH);
}

inline float Vis_GGX(float roughness, float NoV, float NoL)
{
    // Cook and Torrance,
    //  "A Reflectance Model for Computer Graphics"

    const float a = roughness * roughness;
    const float a2 = a * a;

    const float V = NoV + Sqrt( a2 + (NoV - a2 * NoV) * NoV );
    const float L = NoL + Sqrt( a2 + (NoL - a2 * NoL) * NoL );

    return 1.0f / (V * L);
}

inline float Vis_SchlickGGX(float roughness, float NoV, float NoL)
{
    const float a = roughness * roughness;

    const float k = a * 0.5f;

    const float V = NoV - NoV * k + k;
    const float L = NoL - NoL * k + k;

    return 0.25f / (V * L);
}

//------------------------------------------------------------------------------
// Specular F - Fresnel
//------------------------------------------------------------------------------

inline Vec3 F_None(Vec3 const& f0)
{
    return f0;
}

inline Vec3 F_Schlick(Vec3 const& f0, float VoH)
{
    // Schlick,
    //  "An Inexpensive BRDF Model for Physically-Based Rendering"

    const float f = Pow5(1.0f - VoH);

    return f0 - f * f0 + f; // (1.0f - f) * f0 + f; // f0 + (1.0f - f0) * f
}

inline Vec3 F_Schlick_Unreal(Vec3 const& f0, float VoH)
{
    const float f = Pow5(1 - VoH);

    // Anything less than 2% is physically impossible and is instead considered
    // to be shadowing
    return (1.0f - f) * f0 + f * Saturate(50.0f * f0.y);
}

inline Vec3 F_CookTorrance(Vec3 const& f0, float VoH)
{
    const Vec3 s = Sqrt(Clamp(f0, 0.0f, 0.999f));
    const Vec3 n = (1.0f + s) / (1.0f - s);
    const Vec3 g = Sqrt(n * n + VoH * VoH - 1.0f);

    const Vec3 x = (g - VoH) / (g + VoH);
    const Vec3 y = ((g + VoH) * VoH - 1.0f) / ((g - VoH) * VoH + 1);

    return 0.5f * Pow2(x) * (1.0f + Pow2(y));
}

} // namespace brdf
} // namespace math
