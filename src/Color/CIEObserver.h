#pragma once

#define LIB_MATH_COLOR_CIE_OBSERVER_H 1

#include "../Vector.h"

namespace math {

//------------------------------------------------------------------------------
// CIE 1931 2 deg Standard Observer
//------------------------------------------------------------------------------

class CIEObserver
{
    static const Vec3 Table[];

public:
    static Vec3 eval(float lambda)
    {
        if (lambda < 355 || lambda >= 835)
            return { 0, 0, 0 };

        auto i = static_cast<int>(Floor(lambda)) / 5;
        auto s = (lambda - i * 5.0f) / 5.0f;

        i -= 355 / 5;

        return Lerp(Table[i], Table[i + 1], s);
    }

    Vec3 operator()(float lambda) const
    {
        return eval(lambda);
    }
};

//------------------------------------------------------------------------------
// CIE 1964 10 deg Standard Observer
//------------------------------------------------------------------------------

class CIEObserver1964
{
    static const Vec3 Table[];

public:
    static Vec3 eval(float lambda)
    {
        if (lambda < 355 || lambda >= 835)
            return { 0, 0, 0 };

        auto i = static_cast<int>(Floor(lambda)) / 5;
        auto s = (lambda - i * 5.0f) / 5.0f;

        i -= 355 / 5;

        return Lerp(Table[i], Table[i + 1], s);
    }

    Vec3 operator()(float lambda) const
    {
        return eval(lambda);
    }
};

//------------------------------------------------------------------------------
// CIEObserverAnalytic
//
// Color matching functions for CIE 1931 2 deg Standard Observer.
// Using analytic approximations from:
//
// https://research.nvidia.com/publication/simple-analytic-approximations-cie-xyz-color-matching-functions
//------------------------------------------------------------------------------

class CIEObserverAnalytic
{
public:
    static Vec3 eval(float lambda)
    {
        return { CIE_x(lambda), CIE_y(lambda), CIE_z(lambda) };
    }

    Vec3 operator()(float lambda) const
    {
        return eval(lambda);
    }

private:
    static float CIE_x(float lambda)
    {
        float t1 = (lambda - 442.0f) * ((lambda < 442.0f) ? 0.0624f : 0.0374f);
        float t2 = (lambda - 599.8f) * ((lambda < 599.8f) ? 0.0264f : 0.0323f);
        float t3 = (lambda - 501.1f) * ((lambda < 501.1f) ? 0.0490f : 0.0382f);

        return 0.362f * exp(-0.5f * t1 * t1)
             + 1.056f * exp(-0.5f * t2 * t2)
             - 0.065f * exp(-0.5f * t3 * t3);
    }

    static float CIE_y(float lambda)
    {
        float t1 = (lambda - 568.8f) * ((lambda < 568.8f) ? 0.0213f : 0.0247f);
        float t2 = (lambda - 530.9f) * ((lambda < 530.9f) ? 0.0613f : 0.0322f);

        return 0.821f * exp(-0.5f * t1 * t1) + 0.286f * exp(-0.5f * t2 * t2);
    }

    static float CIE_z(float lambda)
    {
        float t1 = (lambda - 437.0f) * ((lambda < 437.0f) ? 0.0845f : 0.0278f);
        float t2 = (lambda - 459.0f) * ((lambda < 459.0f) ? 0.0385f : 0.0725f);

        return 1.217f * exp(-0.5f * t1 * t1) + 0.681f * exp(-0.5f * t2 * t2);
    }
};

} // namespace math
