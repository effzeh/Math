#pragma once

#define LIB_MATH_TMATH_H 1

#include <cstdint>
#include <cmath>

namespace math {

static constexpr double kZero                    =  0.0;
static constexpr double kHalf                    =  0.5;
static constexpr double kOne                     =  1.0;
static constexpr double kNegOne                  = -1.0;
static constexpr double kTwo                     =  2.0;
static constexpr double kSqrtTwo                 =  1.41421356237309504880168872421e+00;
static constexpr double kInvSqrtTwo              =  7.07106781186547524400844362105e-01;
static constexpr double kDegToRad                =  1.74532925199432957692369076849e-02;
static constexpr double kRadToDeg                =  5.72957795130823208767981548141e+01;
static constexpr double kPi                      =  3.14159265358979323846264338328e+00;
static constexpr double kInvPi                   =  3.18309886183790671537767526745e-01;
static constexpr double kPiHalf                  =  1.57079632679489661923132169164e+00;
static constexpr double kTwoPi                   =  6.28318530717958647692528676656e+00;
static constexpr double kFourPi                  =  1.25663706143591729538505735331e+01;
static constexpr double kInvTwoPi                =  1.59154943091895335768883763373e-01;
static constexpr double kInvFourPi               =  7.95774715459476678844418816863e-02;
static constexpr double kSqrtPi                  =  1.77245385090551602729816748334e+00;
static constexpr double kSqrtTwoPi               =  2.50662827463100050241576528481e+00;
static constexpr double kInvSqrtTwoPi            =  3.98942280401432677939946059934e-01;
static constexpr double kFirstRadiationConstant  =  1.191042953e-16;
static constexpr double kSecondRadiationConstant =  1.43877736e-2;

#if 0
template <class T> T Abs(T const& x) { return static_cast<T>(std::abs(x)); }
#else
// Return type is decltype(x < 0 ? -x : x)
inline int                Abs(signed char        x) { return x < 0 ? -x : x; }
inline int                Abs(short              x) { return x < 0 ? -x : x; }
inline int                Abs(int                x) { return x < 0 ? -x : x; } // XXX: PRE: x != MIN
inline long               Abs(long               x) { return x < 0 ? -x : x; } // XXX: PRE: x != MIN
inline long long          Abs(long long          x) { return x < 0 ? -x : x; } // XXX: PRE: x != MIN
inline unsigned int       Abs(unsigned char      x) { return x; }
inline unsigned int       Abs(unsigned short     x) { return x; }
inline unsigned int       Abs(unsigned int       x) { return x; }
inline unsigned long      Abs(unsigned long      x) { return x; }
inline unsigned long long Abs(unsigned long long x) { return x; }
inline float              Abs(float              x) { return std::abs(x); }
inline double             Abs(double             x) { return std::abs(x); }
inline long double        Abs(long double        x) { return std::abs(x); }
#endif

template <class FloatT> FloatT Acos     (FloatT const& x)                  { return std::acos(x);        }
template <class FloatT> FloatT Asin     (FloatT const& x)                  { return std::asin(x);        }
template <class FloatT> FloatT Atan2    (FloatT const& y, FloatT const& x) { return std::atan2(y, x);    }
template <class FloatT> FloatT Atan     (FloatT const& x)                  { return std::atan(x);        }
template <class FloatT> FloatT Ceil     (FloatT const& x)                  { return std::ceil(x);        }
template <class FloatT> FloatT Copysign (FloatT const& x, FloatT const& y) { return std::copysign(x, y); }
template <class FloatT> FloatT Cos      (FloatT const& x)                  { return std::cos(x);         }
template <class FloatT> FloatT Exp      (FloatT const& x)                  { return std::exp(x);         }
template <class FloatT> FloatT Exp2     (FloatT const& x)                  { return std::exp2(x);        }
template <class FloatT> FloatT Floor    (FloatT const& x)                  { return std::floor(x);       }
template <class FloatT> FloatT Fmod     (FloatT const& x, FloatT const& y) { return std::fmod(x, y);     }
template <class FloatT> FloatT Log      (FloatT const& x)                  { return std::log(x);         }
template <class FloatT> FloatT Log2     (FloatT const& x)                  { return std::log2(x);        }
template <class FloatT> FloatT Log10    (FloatT const& x)                  { return std::log10(x);       }
template <class FloatT> FloatT Pow      (FloatT const& x, FloatT const& y) { return std::pow(x, y);      }
template <class FloatT> FloatT Round    (FloatT const& x)                  { return std::round(x);       }
template <class FloatT> FloatT Sin      (FloatT const& x)                  { return std::sin(x);         }
template <class FloatT> FloatT Sqrt     (FloatT const& x)                  { return std::sqrt(x);        }
template <class FloatT> FloatT Tan      (FloatT const& x)                  { return std::tan(x);         }
template <class FloatT> FloatT Trunc    (FloatT const& x)                  { return std::trunc(x);       }

#if 0 // make sure we don't use these...
#define abs ?
#define acos ?
#define asin ?
#define atan2 ?
#define atan ?
#define ceil ?
#define copysign ?
#define cos ?
#define exp ?
#define exp2 ?
#define floor ?
#define fmod ?
#define log ?
#define log2 ?
#define log10 ?
#define pow ?
#define round ?
#define sin ?
#define sqrt ?
#define tan ?
#define trunc ?
#endif

//inline float  Min(float  const& x, float  const& y) { return std::fmin(x, y); }
//inline double Min(double const& x, double const& y) { return std::fmin(x, y); }
//inline float  Max(float  const& x, float  const& y) { return std::fmax(x, y); }
//inline double Max(double const& x, double const& y) { return std::fmax(x, y); }

template <class T>
inline constexpr T const& Min(T const& x, T const& y)
{
    return y < x ? y : x;
}

template <class T>
inline constexpr T const& Max(T const& x, T const& y)
{
    return y < x ? x : y;
}

template <class T, class... Ts>
inline constexpr T const& Min(T const& x, T const& y, Ts const&... z)
{
    return Min(Min(x, y), z...);
}

template <class T, class... Ts>
inline constexpr T const& Max(T const& x, T const& y, Ts const&... z)
{
    return Max(Max(x, y), z...);
}

// Returns Min(x,1)
template <class FloatT>
inline constexpr FloatT Min1(FloatT const& x)
{
    return Min(x, FloatT(1.0));
}

// Returns Max(0,x)
template <class T>
inline constexpr T Max0(T const& x)
{
    return Max(T(0.0), x);
}

// Returns Max(-1,x)
template <class T>
inline constexpr T MaxNeg1(T const& x)
{
    return Max(T(-1.0), x);
}

// Returns x clamped to [lo, hi]
#if 0
template <class T>
inline constexpr T Clamp(T const& x, T const& lo, T const& hi) { return Max(lo, Min(x, hi)); }
#else
template <class T, class U>
inline constexpr T Clamp(T const& x, U const& lo, U const& hi)
{
    return Max(T(lo), Min(x, T(hi)));
}
#endif

// Returns x clamped to [-1, 1]
#if 0
template <class FloatT>
inline constexpr FloatT Clamp(FloatT const& x) { return Clamp(x, FloatT(-1.0), FloatT(1.0)); }
#else
template <class FloatT>
inline constexpr FloatT Clamp(FloatT const& x)
{
    return MaxNeg1(Min1(x));
}
#endif

// Returns x clamped to [0, 1]
template <class FloatT>
inline constexpr FloatT Saturate(FloatT const& x)
{
    return Max0(Min1(x));
}

// Returns the sign of x in {-1,0,1}
template <class T>
inline constexpr T Sign(T const& x)
{
    return x < T(0.0) ? T(-1.0) : (T(0.0) < x ? T(1.0) : T(0.0));
}

// Returns the sign of x in {-1,1}
#if 1
template <class T>
inline constexpr T SignNotZero(T const& x)
{
    return x < T(0.0) ? T(-1.0) : T(1.0);
}
#else
template <class FloatT>
inline FloatT SignNotZero(FloatT const& x)
{
    return Copysign(FloatT(1.0), x);
}
#endif

// Returns 1 / x
template <class FloatT>
inline constexpr FloatT Rcp(FloatT const& x)
{
    return FloatT(1.0) / x;
}

// Returns 1 / sqrt(x)
template <class FloatT>
inline FloatT RcpSqrt(FloatT const& x)
{
    return Sqrt(x) / x;
//  return FloatT(1.0) / Sqrt(x);
}

// Returns x^2
template <class T>
inline constexpr T Square(T const& x)
{
    return x * x;
}

template <class T> inline constexpr T Pow2(T const& x) { return x * x; }
template <class T> inline constexpr T Pow3(T const& x) { return Pow2(x) * x; }
template <class T> inline constexpr T Pow4(T const& x) { return Pow2(Pow2(x)); }
template <class T> inline constexpr T Pow5(T const& x) { return Pow4(x) * x; }

// Returns 1 / x^2
template <class FloatT>
inline constexpr FloatT RcpSquare(FloatT const& x)
{
    return Rcp(Square(x));
}

template <class FloatT>
inline FloatT SignedSqrt(FloatT const& x)
{
    return Copysign(Sqrt(Abs(x)), x);
}

template <class FloatT>
inline FloatT SafeSqrt(FloatT const& x)
{
    return Sqrt(Max0(x));
}

template <class FloatT>
inline FloatT SafeAcos(FloatT const& x)
{
    return Acos(Clamp(x));
}

template <class FloatT>
inline FloatT SafeAsin(FloatT const& x)
{
    return Asin(Clamp(x));
}

// Returns sin(x) / x
template <class FloatT>
FloatT Sinc(FloatT const& x)
{
    return x == FloatT(0.0) ? FloatT(1.0) : Sin(x) / x;
}

template <class FloatT>
void SinCos(FloatT const& phi, FloatT& s, FloatT& c)
{
    s = Sin(phi);
    c = Cos(phi);
}

#ifdef _GNU_SOURCE
inline void SinCos(float const& phi, float& s, float& c)
{
    ::sincosf(phi, &s, &c);
}

inline void SinCos(double const& phi, double& s, double& c)
{
    ::sincos(phi, &s, &c);
}
#endif

// Linear interpolation
template <class T, class FloatT>
inline constexpr T Mix(T const& lo, T const& hi, FloatT const& t)
{
    // Do NOT use: lo + (hi - lo) * t
    // It is inaccurate for t = 1

    return (FloatT(1.0) - t) * lo + t * hi;
}

// Linear interpolation
template <class T, class FloatT>
inline constexpr T Lerp(T const& lo, T const& hi, FloatT const& t)
{
    // Do NOT use: lo + (hi - lo) * t
    // It is inaccurate for t = 1

    return (FloatT(1.0) - t) * lo + t * hi;
}

//template <class T, class U>
//inline constexpr T LerpFast(T const& lo, T const& hi, U const& t) { return lo + t * (hi - lo); }

// Bilinear interpolation
template <class T, class FloatT>
inline constexpr T Bilerp(T const& lo0, T const& hi0, T const& lo1, T const& hi1, FloatT const& t, FloatT const& u)
{
    T u0 = Lerp(lo0, hi0, t);
    T u1 = Lerp(lo1, hi1, t);

    return Lerp(u0, u1, u);
}

template <class FloatT>
inline constexpr FloatT Step(FloatT const& edge, FloatT const& x)
{
    return x < edge ? FloatT(0.0) : FloatT(1.0);
}

template <class FloatT>
inline constexpr FloatT LinearStep(FloatT const& edge0, FloatT const& edge1, FloatT const& x)
{
    return Saturate((x - edge0) / (edge1 - edge0));
}

template <class FloatT>
inline /*constexpr*/ FloatT SmoothStep(FloatT const& edge0, FloatT const& edge1, FloatT const& x)
{
    FloatT y = LinearStep(edge0, edge1, x);
    return y * y * (FloatT(3.0) - FloatT(2.0) * y);
}

// Integer square root.
template <class IntT>
inline /*constexpr*/ IntT IntSqrt(IntT n)
{
    if (n <= IntT(0))
        return IntT(0);

    IntT m = n;
    IntT k = IntT(1);
    do
    {
        m = k + (m - k) / IntT(2);
        k = n / m;
    } while (m > k);

    return m;
}

// Returns the determinant of [m00 m01; m10 m11]
//
// NOTE:
// Row-major!
//
template <class T>
inline constexpr T Det2(T const& m00, T const& m01, T const& m10, T const& m11)
{
    return m00 * m11 - m10 * m01;
}

// Returns real roots for b x + c in t.
// Returns the number of real roots.
template <class FloatT>
int SolveLinear(FloatT const& b, FloatT const& c, FloatT& t)
{
    if (b == FloatT(0.0))
        return 0;

    t = -c / b;

    return 1;
}

// Returns real roots for a x^2 + b x + c in t1 and t2.
// Returns the number of real roots.
template <class FloatT>
int SolveQuadratic(FloatT const& a, FloatT const& b, FloatT const& c, FloatT& t1, FloatT& t2)
{
    if (a == FloatT(0.0))
        return SolveLinear(b, c, t1);

    auto discr = b * b - FloatT(4.0) * a * c;
    if (discr > FloatT(0.0))
    {
        auto t = FloatT(-0.5) * (b + (b >= FloatT(0.0) ? Sqrt(discr) : -Sqrt(discr)));

        t1 = t / a;
        t2 = c / t;

        return 2;
    }

    if (discr < FloatT(0.0))
        return 0;

    t1 = FloatT(-0.5) * b / a;

    return 1;
}

#if 1

// Returns: condition < 0 ? x : y
inline constexpr float SelectNegative(float condition, float x, float y)
{
#if 1
    return (condition < 0) ? x : y;
#else
    //
    // result
    //  = condition < 0 ? x : y
    //  = condition < 0 ? ((1 - 0) * x + 0 * y) : ((1 - 1) * x + 1 * y)
    //  = condition < 0 ? Mix(x, y, 0) : Mix(x, y, 1)
    //  = Mix(x, y, Step(0, condition))
    //
    return Mix(x, y, Step(0, condition));
#endif
}

// Returns: condition > 0 ? x : y
inline constexpr float SelectPositive(float condition, float x, float y)
{
#if 1
    return (0 < condition) ? x : y;
#else
    //
    // result
    //  = condition > 0 ? x : y
    //  = condition > 0 ? ((1 - 0) * x + 0 * y) : ((1 - 1) * x + 1 * y)
    //  = condition > 0 ? Mix(x, y, 0) : Mix(x, y, 1)
    //  = Mix(x, y, Step(condition, 0))
    //
    return Mix(x, y, Step(condition, 0));
#endif
}

#endif

#if 0

inline float IntAsFloat(int n)
{
    union { int n; float f; } bits = {n};
    return bits.f;
}

inline int FloatAsInt(float f)
{
    union { float f; int n; } bits = {f};
    return bits.n;
}

// Returns the signbit of X encoded as +-0
inline float SignMask(float x) { return IntAsFloat( FloatAsInt(x) & 0x80000000 ); }

inline float BitAnd(float x, float y) { return IntAsFloat( FloatAsInt(x) & FloatAsInt(y) ); }
inline float BitOr (float x, float y) { return IntAsFloat( FloatAsInt(x) | FloatAsInt(y) ); }
inline float BitXor(float x, float y) { return IntAsFloat( FloatAsInt(x) ^ FloatAsInt(y) ); }

#endif

// Returns x/y rounded up
template <class IntT>
inline constexpr IntT DivUp(IntT x, IntT y)
{
    return (x + (y - IntT(1))) / y;
}

// Round up X to the nearest multiple of Y
template <class IntT>
inline constexpr IntT RoundUp(IntT const& x, IntT const& y)
{
    return DivUp(x, y) * y;
}

// Aligns N up to the nearest multiple of A.
// A must be a power of 2.
template <class IntT>
inline constexpr IntT AlignUp(IntT n, IntT a)
{
    return (n + (a - IntT(1))) & ~(a - IntT(1));
}

// Return true if x is in [a, b)
template <class T>
inline constexpr bool InRange(T const& x, T const& a, T const& b)
{
    //  return a <= x && x < b;
    return !(x < a) && x < b;
}

// Return true if x is in [a, b]
template <class T>
inline constexpr bool InRangeInclusive(T const& x, T const& a, T const& b)
{
    //  return a <= x && x <= b;
    //  return !(x < a) && !(b < x);
    return !((x < a) || (b < x));
}

#if 0

// Normalize angle into the range 0 <= X < 2pi.
template <class FloatT>
FloatT NormalizeAngle2Pi(FloatT x)
{
    FloatT w = Fmod(x, FloatT(kTwoPi));

    if (w < FloatT(0.0))
        w += FloatT(kTwoPi);

    return w;
}

// Normalize angle into the range -pi <= X < pi.
template <class FloatT>
FloatT NormalizeAnglePi(FloatT x)
{
    FloatT w = Fmod(x, FloatT(kTwoPi));

    if (Abs(w) >= FloatT(kPi))
        w -= Copysign(FloatT(kPi), x);

    return w;
}

#endif

// Computes sqrt(x^2 + y^2) avoiding unneccesary overflow.
template <class FloatT>
FloatT Hypot(FloatT x, FloatT y)
{
    // LAPACK routine slapy2

    FloatT ax = Abs(x);
    FloatT ay = Abs(y);

#if 1
    // Compute Max/Min the FORTRAN way...
    FloatT w = ax >= ay ? ax : ay;
    FloatT z = ax <= ay ? ax : ay;
#else
    FloatT w = Max(ax, ay);
    FloatT z = Min(ax, ay);
#endif
    if (z == FloatT(0.0))
    {
        return w;
    }

    return w * Sqrt(FloatT(1.0) + Square(z / w));
}

// Returns sqrt(x^2 + y^2 + z^2), taking care not to cause unnecessary overflow.
template <class FloatT>
FloatT Hypot(FloatT x, FloatT y, FloatT z)
{
    // LAPACK routine slapy3

    FloatT ax = Abs(x);
    FloatT ay = Abs(y);
    FloatT az = Abs(z);

#if 1
    // Compute Max/Min the FORTRAN way...
    FloatT s = ax;
    s = (s >= ay) ? s : ay;
    s = (s >= az) ? s : az;
#else
    FloatT s = Max(ax, ay, az);
#endif
    if (s == FloatT(0.0))
    {
        // W can be zero for Max(0,nan,0).
        // Adding all three entries together will make sure NaN will not disappear.
        return ax + ay + az;
    }

    return s * Sqrt(Square(ax / s) + Square(ay / s) + Square(az / s));
}

// Returns sqrt(x^2 + y^2 + z^2 + w^2), taking care not to cause unnecessary overflow.
template <class FloatT>
FloatT Hypot(FloatT x, FloatT y, FloatT z, FloatT w)
{
    FloatT ax = Abs(x);
    FloatT ay = Abs(y);
    FloatT az = Abs(z);
    FloatT aw = Abs(w);

#if 1
    // Compute Max/Min the FORTRAN way...
    FloatT s = ax;
    s = (s >= ay) ? s : ay;
    s = (s >= az) ? s : az;
    s = (s >= aw) ? s : aw;
#else
    FloatT s = Max(ax, ay, az, aw);
#endif
    if (s == FloatT(0.0))
    {
        // W can be zero for Max(0,nan,0).
        // Adding all three entries together will make sure NaN will not disappear.
        return ax + ay + az + aw;
    }

    return s * Sqrt(Square(ax / s) + Square(ay / s) + Square(az / s) + Square(aw / s));
}

// Returns sqrt(x' * x), takingn care not to cause unncessary overflow.
template <class FloatT>
FloatT Norm2(FloatT const* x, int n, int incx = 1)
{
    // LAPACK routine snrm2

    const FloatT zero(0.0);
    const FloatT one(1.0);

    if (n <= 0 || incx <= 0)
        return zero;

    if (n == 1)
        return Abs(x[0]);

    FloatT scale = zero;
    FloatT ssq = one;

    for (int i = 0, iend = n * incx; i < iend; i += incx)
    {
        FloatT a = x[i];

        if (a == zero)
            continue;

        a = Abs(a);
        if (scale < a)
        {
            ssq = one + ssq * Square(scale / a);
            scale = a;
        }
        else
        {
            ssq += Square(a / scale);
        }
    }

    return scale * Sqrt(ssq);
}

inline constexpr int _inc(int i) { return ++i; }
inline constexpr int _inc_nz(int i) { return i + (i != 0); }

// Returns the index of the min-element
template <class T>
inline constexpr int MinIndex(T const& x, T const& y)
{
    return y < x ? 1 : 0;
}

// Returns the index of the min-element
template <class T, class... Ts>
inline constexpr int MinIndex(T const& x, T const& y, Ts const&... ints)
{
    return y < x
        ? 1 + MinIndex(y, ints...)
        : _inc_nz(MinIndex(x, ints...));
}

// Returns the index of the max-element
template <class T>
inline constexpr int MaxIndex(T const& x, T const& y)
{
    return y < x ? 0 : 1;
}

// Returns the index of the max-element
template <class T, class... Ts>
inline constexpr int MaxIndex(T const& x, T const& y, Ts const&... ints)
{
    return y < x
        ? _inc_nz(MaxIndex(x, ints...))
        : 1 + MaxIndex(y, ints...);
}

// Mask out the bits which are 0 in prev and 1 in curr
template <class UnsignedIntT>
inline constexpr UnsignedIntT RecentlySet(UnsignedIntT prev, UnsignedIntT curr)
{
    return curr & ~(prev & curr);
}

// Mask out the bits which are 1 in prev and 0 in curr
template <class UnsignedIntT>
inline constexpr UnsignedIntT RecentlyCleared(UnsignedIntT prev, UnsignedIntT curr)
{
    return prev & ~(prev & curr);
}

} // namespace math
