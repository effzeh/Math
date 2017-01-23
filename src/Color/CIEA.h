#pragma once

#define LIB_MATH_COLOR_CIE_A_H 1

#include <cmath>

namespace math {

// Spectral Power Distribution of CIE A-Illuminant.
// Normalized such that P(560) = 1.
class CIEA
{
public:
    // See:
    // https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_A

    float S = 8082.192095566608524471400990977978939540796622352193037455f;

public:
    float operator()(float lambda) const
    {
        return ::powf(560.0f / lambda, 5.0f) * (S / eval(lambda));
    }

private:
    static float eval(float lambda)
    {
        return ::expf(0.01435e+9f / (2848.0f * lambda)) - 1.0f;
    }
};

} // namespace math
