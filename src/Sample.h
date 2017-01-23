#pragma once

#define LIB_MATH_SAMPLE_H 1

#include "Vector.h"

namespace math {

struct Sample3
{
    // Location of the sample
    Vec3 value;
    // Probability density of the sample
    float pdf;

    Sample3() = default;
    Sample3(Vec3 const& value, float pdf) : value(value), pdf(pdf) {}

    static Sample3 SphereUniform(float u, float v)
    {
        const float phi       = u * float{kTwoPi};
        const float cos_theta = 2.0f * v - 1.0f;

        return { Vec3::Spherical(phi, cos_theta), float{kInvFourPi} };
    }

    static Sample3 SphereCosine(float u, float v)
    {
        const float phi       = u * float{kTwoPi};
        const float cos_theta = Sqrt(Abs(2.0f * v - 1.0f));

        return { Vec3::Spherical(phi, cos_theta), 2.0f * float{kInvPi} * cos_theta };
    }

    static Sample3 Hemisphere(float u, float v, float exp, float angle)
    {
        // assert: exp >= 0
        // assert: angle > 0
        // assert: angle <= pi/2

        const float K         = 1.0f - Pow(Cos(angle), exp + 1.0f);

        const float phi       = u * float{kTwoPi};
        const float cos_theta = Pow(Saturate(1.0f - v * K), 1.0f / (exp + 1.0f));
        const float pdf       = float{kInvTwoPi} * (exp + 1.0f) * Pow(cos_theta, exp) / K;

        return { Vec3::Spherical(phi, cos_theta), pdf };
    }

    static Sample3 HemisphereUniform(float u, float v)
    {
        // exp = 0, angle = pi/2
        //
        // K = 1
        // and use v ~ 1-v

        const float phi       = u * float{kTwoPi};
        const float cos_theta = v;

        return { Vec3::Spherical(phi, cos_theta), float{kInvTwoPi} };
    }

    static Sample3 HemisphereCosine(float u, float v)
    {
        // exp = 1, angle = pi/2
        //
        // K = 1
        // and use v ~ 1-v

        const float phi       = u * float{kTwoPi};
        const float cos_theta = Sqrt(v);
        const float pdf       = float{kInvPi} * cos_theta;

        return { Vec3::Spherical(phi, cos_theta), pdf };
    }

    static Sample3 HemisphereCosinePower(float u, float v, float exp)
    {
        // angle = pi/2
        //
        // K = 1
        // and use v ~ 1-v

        const float phi       = u * float{kTwoPi};
        const float cos_theta = Pow(v, 1.0f / (exp + 1.0f));
        const float pdf       = float{kInvTwoPi} * (exp + 1.0f) * Pow(cos_theta, exp);

        return { Vec3::Spherical(phi, cos_theta), pdf };
    }

    static Sample3 ConeUniform(float u, float v, float angle)
    {
        // exp = 0

        const float K         = 1.0f - Cos(angle);

        const float phi       = u * float{kTwoPi};
        const float cos_theta = Saturate(1.0f - v * K);
        const float pdf       = float{kInvTwoPi} / K;

        return { Vec3::Spherical(phi, cos_theta), pdf };
    }

    static Sample3 ConeCosine(float u, float v, float angle)
    {
        // exp = 1

        const float K         = 1.0f - Square(Cos(angle));

        const float phi       = u * float{kTwoPi};
        const float cos_theta = Sqrt(Saturate(1.0f - v * K));
        const float pdf       = cos_theta * float{kInvPi} / K;

        return { Vec3::Spherical(phi, cos_theta), pdf };
    }
};

} // namespace math
