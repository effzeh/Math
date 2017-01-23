#pragma once

#define LIB_MATH_COLOR_BLACKBODY_H 1

#include <cmath>

namespace math {

// Spectral Power Distribution of a black body
// Normalized such that P(560) = 1.
class Blackbody
{
    // See:
    // https://en.wikipedia.org/wiki/Planckian_locus

    float T; // Temperature (Kelvin)
    float S; // Normalization factor

public:
    Blackbody() : Blackbody(6504.0f)
    {
    }

    explicit Blackbody(float temp)
        : T(temp)
        , S(eval(T, 560.0f))
    {
    }

    float operator()(float lambda) const
    {
        return ::powf(560.0f / lambda, 5.0f) * (S / eval(T, lambda));
    }

private:
    static float eval(float T, float lambda)
    {
        // Second radiation constant (from CODATA 2014)
        //  c2 = 0.0143877736
        // The 1/10^-9 is to convert nm to m.
        return ::expf(0.0143877736e+9f / (T * lambda)) - 1.0f;
    }
};

} // namespace math
