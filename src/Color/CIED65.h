#pragma once

#define LIB_MATH_COLOR_CIE_D65_H 1

#include <cmath>

namespace math {

// Spectral Power Distribution of CIE D65-Illuminant.
// Normalized such that P(560) = 1.
class CIED65
{
    static const float Table[];

public:
    static float eval(float lambda)
    {
        if (lambda < 295 || lambda >= 835)
            return 0;

        auto i = static_cast<int>(::floorf(lambda)) / 5;
        auto s = (lambda - i * 5.0f) / 5.0f;

        i -= 295 / 5;

        return (1.0f - s) * Table[i] + s * Table[i + 1];
    }

    float operator()(float lambda) const
    {
        return eval(lambda);
    }
};

} // namespace math
