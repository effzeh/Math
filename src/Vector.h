#pragma once

#define LIB_MATH_VECTOR_H 1

#include "TMath.h"

namespace math {

//==============================================================================
// Class definitions
//==============================================================================

struct Vec2b;
struct Vec2;
struct Vec2i;
struct Vec3b;
struct Vec3;
struct Vec3i;
struct Vec4b;
struct Vec4;
struct Vec4i;
struct Point3; // (x,y,z,1)

struct Vec2b
{
    bool x;
    bool y;

    Vec2b() = default;

    Vec2b(bool x, bool y) : x(x), y(y) {}

    explicit Vec2b(bool t) : x(t), y(t) {}

    bool& operator[](int i) { return (&x)[i]; }

    bool const& operator[](int i) const { return (&x)[i]; }
};

struct Vec2
{
    float x;
    float y;

    Vec2() = default;

    Vec2(float x, float y) : x(x), y(y) {}

    explicit Vec2(float t) : x(t), y(t) {}

    explicit Vec2(float const* v) : x(v[0]), y(v[1]) {}

    explicit Vec2(Vec2i const& v);
    explicit Vec2(Vec3 const& v);
    explicit Vec2(Vec3i const& v);
    explicit Vec2(Vec4 const& v);
    explicit Vec2(Vec4i const& v);

    float* data() { return &x; }

    float const* data() const { return &x; }

    float& operator[](int i) { return (&x)[i]; }

    float const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Vec2 FromCoords(V const& v);
};

struct Vec2i
{
    int x;
    int y;

    Vec2i() = default;

    Vec2i(int x, int y) : x(x), y(y) {}

    explicit Vec2i(int t) : x(t), y(t) {}

    explicit Vec2i(int const* v) : x(v[0]), y(v[1]) {}

    explicit Vec2i(Vec2 const& v);
    explicit Vec2i(Vec3 const& v);
    explicit Vec2i(Vec3i const& v);
    explicit Vec2i(Vec4 const& v);
    explicit Vec2i(Vec4i const& v);

    int* data() { return &x; }

    int const* data() const { return &x; }

    int& operator[](int i) { return (&x)[i]; }

    int const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Vec2i FromCoords(V const& v);
};

struct Vec3b
{
    bool x;
    bool y;
    bool z;

    Vec3b() = default;

    Vec3b(bool x, bool y, bool z) : x(x), y(y), z(z) {}

    explicit Vec3b(bool t) : x(t), y(t), z(t) {}

    bool& operator[](int i) { return (&x)[i]; }

    bool const& operator[](int i) const { return (&x)[i]; }
};

struct Vec3
{
    float x;
    float y;
    float z;

    Vec3() = default;

    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    explicit Vec3(float t) : x(t), y(t), z(t) {}

    explicit Vec3(float const* v) : x(v[0]), y(v[1]), z(v[2]) {}

    explicit Vec3(Vec2 const& v, float z = 0);
    explicit Vec3(Vec2i const& v, float z = 0);
    explicit Vec3(Vec3i const& v);
    explicit Vec3(Vec4 const& v);
    explicit Vec3(Vec4i const& v);
    explicit Vec3(Point3 const& p);

    float* data() { return &x; }

    float const* data() const { return &x; }

    float& operator[](int i) { return (&x)[i]; }

    float const& operator[](int i) const { return (&x)[i]; }

    // Construct from spherical coordinates
    static Vec3 Spherical(float phi, float cosTheta, float sinTheta);

    // Construct from spherical coordinates
    static Vec3 Spherical(float phi, float cosTheta);

    template <class V>
    static Vec3 FromCoords(V const& v);
};

struct Vec3i
{
    int x;
    int y;
    int z;

    Vec3i() = default;

    Vec3i(int x, int y, int z) : x(x), y(y), z(z) {}

    explicit Vec3i(int t) : x(t), y(t), z(t) {}

    explicit Vec3i(int const* v) : x(v[0]), y(v[1]), z(v[2]) {}

    explicit Vec3i(Vec2 const& v, int z = 0);
    explicit Vec3i(Vec2i const& v, int z = 0);
    explicit Vec3i(Vec3 const& v);
    explicit Vec3i(Vec4 const& v);
    explicit Vec3i(Vec4i const& v);

    int* data() { return &x; }

    int const* data() const { return &x; }

    int& operator[](int i) { return (&x)[i]; }

    int const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Vec3i FromCoords(V const& v);
};

struct Vec4b
{
    bool x;
    bool y;
    bool z;
    bool w;

    Vec4b() = default;

    Vec4b(bool x, bool y, bool z, bool w) : x(x), y(y), z(z), w(w) {}

    explicit Vec4b(bool t) : x(t), y(t), z(t), w(t) {}

    bool& operator[](int i) { return (&x)[i]; }

    bool const& operator[](int i) const { return (&x)[i]; }
};

struct Vec4
{
    float x;
    float y;
    float z;
    float w;

    Vec4() = default;

    Vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

    explicit Vec4(float t) : x(t), y(t), z(t), w(t) {}

    explicit Vec4(float const* v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}

    explicit Vec4(Vec2 const& v, float z = 0, float w = 0);
    explicit Vec4(Vec2i const& v, float z = 0, float w = 0);
    explicit Vec4(Vec3 const& v, float w = 0);
    explicit Vec4(Vec3i const& v, float w = 0);
    explicit Vec4(Vec4i const& v);
    explicit Vec4(Point3 const& p);

    float* data() { return &x; }

    float const* data() const { return &x; }

    float& operator[](int i) { return (&x)[i]; }

    float const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Vec4 FromCoords(V const& v);
};

struct Vec4i
{
    int x;
    int y;
    int z;
    int w;

    Vec4i() = default;

    Vec4i(int x, int y, int z, int w) : x(x), y(y), z(z), w(w) {}

    explicit Vec4i(int t) : x(t), y(t), z(t), w(t) {}

    explicit Vec4i(int const* v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}

    explicit Vec4i(Vec2 const& v, int z = 0, int w = 0);
    explicit Vec4i(Vec2i const& v, int z = 0, int w = 0);
    explicit Vec4i(Vec3 const& v, int w = 0);
    explicit Vec4i(Vec3i const& v, int w = 0);
    explicit Vec4i(Vec4 const& v);

    int* data() { return &x; }

    int const* data() const { return &x; }

    int& operator[](int i) { return (&x)[i]; }

    int const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Vec4i FromCoords(V const& v);
};

struct Point3
{
    float x;
    float y;
    float z;

    Point3() = default;

    Point3(float x, float y, float z) : x(x), y(y), z(z) {}

    explicit Point3(float t) : x(t), y(t), z(t) {}

    explicit Point3(float const* v) : x(v[0]), y(v[1]), z(v[2]) {}

    // Construct from 3D vector
    explicit Point3(Vec3 const& v);
    // Construct from 4D vector (performs perspective division)
    explicit Point3(Vec4 const& v);

    float* data() { return &x; }

    float const* data() const { return &x; }

    float& operator[](int i) { return (&x)[i]; }

    float const& operator[](int i) const { return (&x)[i]; }

    template <class V>
    static Point3 FromCoords(V const& v);
};

//==============================================================================
// Vec2b
//==============================================================================

inline bool All (Vec2b const& v) { return v.x && v.y; }
inline bool Any (Vec2b const& v) { return v.x || v.y; }
inline bool None(Vec2b const& v) { return !Any(v); }

inline Vec2b operator&(Vec2b const& u, Vec2b const& v) { return { u.x && v.x, u.y && v.y }; }
inline Vec2b operator|(Vec2b const& u, Vec2b const& v) { return { u.x || v.x, u.y || v.y }; }

//==============================================================================
// Vec2
//==============================================================================

inline Vec2::Vec2(Vec2i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
{
}

inline Vec2::Vec2(Vec3 const& v) : x(v.x), y(v.y)
{
}

inline Vec2::Vec2(Vec3i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
{
}

inline Vec2::Vec2(Vec4 const& v) : x(v.x), y(v.y)
{
}

inline Vec2::Vec2(Vec4i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
{
}

template <class V>
inline Vec2 Vec2::FromCoords(V const& v) {
    return { static_cast<float>(v.x), static_cast<float>(v.y) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec2 const& u, Vec2 const& v) {
    return u.x == v.x && u.y == v.y;
}

inline bool operator!=(Vec2 const& u, Vec2 const& v) {
    return !(u == v);
}

inline bool operator<(Vec2 const& u, Vec2 const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;

    return u.y < v.y;
}

inline bool operator>(Vec2 const& u, Vec2 const& v) {
    return v < u;
}

inline bool operator<=(Vec2 const& u, Vec2 const& v) {
    return !(v < u);
}

inline bool operator>=(Vec2 const& u, Vec2 const& v) {
    return !(u < v);
}

inline Vec2b CmpEQ(Vec2 const& u, Vec2 const& v) { return { u.x == v.x, u.y == v.y }; }
inline Vec2b CmpNE(Vec2 const& u, Vec2 const& v) { return { u.x != v.x, u.y != v.y }; }
inline Vec2b CmpLT(Vec2 const& u, Vec2 const& v) { return { u.x <  v.x, u.y <  v.y }; }
inline Vec2b CmpLE(Vec2 const& u, Vec2 const& v) { return { u.x <= v.x, u.y <= v.y }; }
inline Vec2b CmpGT(Vec2 const& u, Vec2 const& v) { return { u.x >  v.x, u.y >  v.y }; }
inline Vec2b CmpGE(Vec2 const& u, Vec2 const& v) { return { u.x >= v.x, u.y >= v.y }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec2 operator+(Vec2 const& u) { return u; }
inline Vec2 operator-(Vec2 const& u) { return { -u.x, -u.y }; }

inline Vec2 operator+(Vec2 const& u, Vec2 const& v) { return { u.x + v.x, u.y + v.y }; }
inline Vec2 operator-(Vec2 const& u, Vec2 const& v) { return { u.x - v.x, u.y - v.y }; }
inline Vec2 operator*(Vec2 const& u, Vec2 const& v) { return { u.x * v.x, u.y * v.y }; }
inline Vec2 operator/(Vec2 const& u, Vec2 const& v) { return { u.x / v.x, u.y / v.y }; }

inline Vec2& operator+=(Vec2& u, Vec2 const& v) { u = u + v; return u; }
inline Vec2& operator-=(Vec2& u, Vec2 const& v) { u = u - v; return u; }
inline Vec2& operator*=(Vec2& u, Vec2 const& v) { u = u * v; return u; }
inline Vec2& operator/=(Vec2& u, Vec2 const& v) { u = u / v; return u; }

inline Vec2 operator+(float u, Vec2 const& v) { return { u + v.x, u + v.y }; }
inline Vec2 operator-(float u, Vec2 const& v) { return { u - v.x, u - v.y }; }
inline Vec2 operator*(float u, Vec2 const& v) { return { u * v.x, u * v.y }; }
inline Vec2 operator/(float u, Vec2 const& v) { return { u / v.x, u / v.y }; }

inline Vec2 operator+(Vec2 const& u, float v) { return { u.x + v, u.y + v }; }
inline Vec2 operator-(Vec2 const& u, float v) { return { u.x - v, u.y - v }; }
inline Vec2 operator*(Vec2 const& u, float v) { return { u.x * v, u.y * v }; }
inline Vec2 operator/(Vec2 const& u, float v) { return { u.x / v, u.y / v }; }

inline Vec2& operator+=(Vec2& u, float v) { u = u + v; return u; }
inline Vec2& operator-=(Vec2& u, float v) { u = u - v; return u; }
inline Vec2& operator*=(Vec2& u, float v) { u = u * v; return u; }
inline Vec2& operator/=(Vec2& u, float v) { u = u / v; return u; }

//------------------------------------------------------------------------------
// Conditional
//

inline Vec2 Select(Vec2b const& condition, Vec2 const& true_, Vec2 const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y };
}

//==============================================================================
// Vec2i
//==============================================================================

inline Vec2i::Vec2i(Vec2 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
{
}

inline Vec2i::Vec2i(Vec3 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
{
}

inline Vec2i::Vec2i(Vec3i const& v)
    : x(v.x)
    , y(v.y)
{
}

inline Vec2i::Vec2i(Vec4 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
{
}

inline Vec2i::Vec2i(Vec4i const& v)
    : x(v.x)
    , y(v.y)
{
}

template <class V>
inline Vec2i Vec2i::FromCoords(V const& v) {
    return { static_cast<int>(v.x), static_cast<int>(v.y) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec2i const& u, Vec2i const& v) {
    return u.x == v.x && u.y == v.y;
}

inline bool operator!=(Vec2i const& u, Vec2i const& v) {
    return !(u == v);
}

inline bool operator<(Vec2i const& u, Vec2i const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;

    return u.y < v.y;
}

inline bool operator>(Vec2i const& u, Vec2i const& v) {
    return v < u;
}

inline bool operator<=(Vec2i const& u, Vec2i const& v) {
    return !(v < u);
}

inline bool operator>=(Vec2i const& u, Vec2i const& v) {
    return !(u < v);
}

inline Vec2b CmpEQ(Vec2i const& u, Vec2i const& v) { return { u.x == v.x, u.y == v.y }; }
inline Vec2b CmpNE(Vec2i const& u, Vec2i const& v) { return { u.x != v.x, u.y != v.y }; }
inline Vec2b CmpLT(Vec2i const& u, Vec2i const& v) { return { u.x <  v.x, u.y <  v.y }; }
inline Vec2b CmpLE(Vec2i const& u, Vec2i const& v) { return { u.x <= v.x, u.y <= v.y }; }
inline Vec2b CmpGT(Vec2i const& u, Vec2i const& v) { return { u.x >  v.x, u.y >  v.y }; }
inline Vec2b CmpGE(Vec2i const& u, Vec2i const& v) { return { u.x >= v.x, u.y >= v.y }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec2i operator+(Vec2i const& u) { return u; }
inline Vec2i operator-(Vec2i const& u) { return { -u.x, -u.y }; }
inline Vec2i operator!(Vec2i const& u) { return { !u.x, !u.y }; }
inline Vec2i operator~(Vec2i const& u) { return { ~u.x, ~u.y }; }

inline Vec2i operator+ (Vec2i const& u, Vec2i const& v) { return { u.x +  v.x, u.y +  v.y }; }
inline Vec2i operator- (Vec2i const& u, Vec2i const& v) { return { u.x -  v.x, u.y -  v.y }; }
inline Vec2i operator* (Vec2i const& u, Vec2i const& v) { return { u.x *  v.x, u.y *  v.y }; }
inline Vec2i operator/ (Vec2i const& u, Vec2i const& v) { return { u.x /  v.x, u.y /  v.y }; }
inline Vec2i operator% (Vec2i const& u, Vec2i const& v) { return { u.x %  v.x, u.y %  v.y }; }
inline Vec2i operator& (Vec2i const& u, Vec2i const& v) { return { u.x &  v.x, u.y &  v.y }; }
inline Vec2i operator| (Vec2i const& u, Vec2i const& v) { return { u.x |  v.x, u.y |  v.y }; }
inline Vec2i operator^ (Vec2i const& u, Vec2i const& v) { return { u.x ^  v.x, u.y ^  v.y }; }
inline Vec2i operator<<(Vec2i const& u, Vec2i const& v) { return { u.x << v.x, u.y << v.y }; }
inline Vec2i operator>>(Vec2i const& u, Vec2i const& v) { return { u.x >> v.x, u.y >> v.y }; }

inline Vec2i& operator+= (Vec2i& u, Vec2i const& v) { u = u +  v; return u; }
inline Vec2i& operator-= (Vec2i& u, Vec2i const& v) { u = u -  v; return u; }
inline Vec2i& operator*= (Vec2i& u, Vec2i const& v) { u = u *  v; return u; }
inline Vec2i& operator/= (Vec2i& u, Vec2i const& v) { u = u /  v; return u; }
inline Vec2i& operator%= (Vec2i& u, Vec2i const& v) { u = u %  v; return u; }
inline Vec2i& operator&= (Vec2i& u, Vec2i const& v) { u = u &  v; return u; }
inline Vec2i& operator|= (Vec2i& u, Vec2i const& v) { u = u |  v; return u; }
inline Vec2i& operator^= (Vec2i& u, Vec2i const& v) { u = u ^  v; return u; }
inline Vec2i& operator<<=(Vec2i& u, Vec2i const& v) { u = u << v; return u; }
inline Vec2i& operator>>=(Vec2i& u, Vec2i const& v) { u = u >> v; return u; }

inline Vec2i operator+ (int u, Vec2i const& v) { return { u +  v.x, u +  v.y }; }
inline Vec2i operator- (int u, Vec2i const& v) { return { u -  v.x, u -  v.y }; }
inline Vec2i operator* (int u, Vec2i const& v) { return { u *  v.x, u *  v.y }; }
inline Vec2i operator/ (int u, Vec2i const& v) { return { u /  v.x, u /  v.y }; }
inline Vec2i operator% (int u, Vec2i const& v) { return { u %  v.x, u %  v.y }; }
inline Vec2i operator& (int u, Vec2i const& v) { return { u &  v.x, u &  v.y }; }
inline Vec2i operator| (int u, Vec2i const& v) { return { u |  v.x, u |  v.y }; }
inline Vec2i operator^ (int u, Vec2i const& v) { return { u ^  v.x, u ^  v.y }; }
inline Vec2i operator<<(int u, Vec2i const& v) { return { u << v.x, u << v.y }; }
inline Vec2i operator>>(int u, Vec2i const& v) { return { u >> v.x, u >> v.y }; }

inline Vec2i operator+ (Vec2i const& u, int v) { return { u.x +  v, u.y +  v }; }
inline Vec2i operator- (Vec2i const& u, int v) { return { u.x -  v, u.y -  v }; }
inline Vec2i operator* (Vec2i const& u, int v) { return { u.x *  v, u.y *  v }; }
inline Vec2i operator/ (Vec2i const& u, int v) { return { u.x /  v, u.y /  v }; }
inline Vec2i operator% (Vec2i const& u, int v) { return { u.x %  v, u.y %  v }; }
inline Vec2i operator& (Vec2i const& u, int v) { return { u.x &  v, u.y &  v }; }
inline Vec2i operator| (Vec2i const& u, int v) { return { u.x |  v, u.y |  v }; }
inline Vec2i operator^ (Vec2i const& u, int v) { return { u.x ^  v, u.y ^  v }; }
inline Vec2i operator<<(Vec2i const& u, int v) { return { u.x << v, u.y << v }; }
inline Vec2i operator>>(Vec2i const& u, int v) { return { u.x >> v, u.y >> v }; }

inline Vec2i& operator+= (Vec2i& u, int v) { u = u +  v; return u; }
inline Vec2i& operator-= (Vec2i& u, int v) { u = u -  v; return u; }
inline Vec2i& operator*= (Vec2i& u, int v) { u = u *  v; return u; }
inline Vec2i& operator/= (Vec2i& u, int v) { u = u /  v; return u; }
inline Vec2i& operator%= (Vec2i& u, int v) { u = u %  v; return u; }
inline Vec2i& operator&= (Vec2i& u, int v) { u = u &  v; return u; }
inline Vec2i& operator|= (Vec2i& u, int v) { u = u |  v; return u; }
inline Vec2i& operator^= (Vec2i& u, int v) { u = u ^  v; return u; }
inline Vec2i& operator<<=(Vec2i& u, int v) { u = u << v; return u; }
inline Vec2i& operator>>=(Vec2i& u, int v) { u = u >> v; return u; }

//------------------------------------------------------------------------------
// Conditional
//

inline Vec2i Select(Vec2b const& condition, Vec2i const& true_, Vec2i const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y };
}

//==============================================================================
// Vec3b
//==============================================================================

inline bool All (Vec3b const& v) { return v.x && v.y && v.z; }
inline bool Any (Vec3b const& v) { return v.x || v.y || v.z; }
inline bool None(Vec3b const& v) { return !Any(v); }

inline Vec3b operator&(Vec3b const& u, Vec3b const& v) { return { u.x && v.x, u.y && v.y, u.z && v.z }; }
inline Vec3b operator|(Vec3b const& u, Vec3b const& v) { return { u.x || v.x, u.y || v.y, u.z || v.z }; }

//==============================================================================
// Vec3
//==============================================================================

inline Vec3::Vec3(Vec2 const& v, float z)
    : x(v.x)
    , y(v.y)
    , z(z)
{
}

inline Vec3::Vec3(Vec2i const& v, float z)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(z)
{
}

inline Vec3::Vec3(Vec3i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(static_cast<float>(v.z))
{
}

inline Vec3::Vec3(Vec4 const& v)
    : x(v.x)
    , y(v.y)
    , z(v.z)
{
}

inline Vec3::Vec3(Vec4i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(static_cast<float>(v.z))
{
}

inline Vec3::Vec3(Point3 const& p)
    : x(p.x)
    , y(p.y)
    , z(p.z)
{
}

inline Vec3 Vec3::Spherical(float phi, float cosTheta, float sinTheta)
{
    auto x = sinTheta * Cos(phi);
    auto y = sinTheta * Sin(phi);
    auto z = cosTheta;

    return {x, y, z};
}

inline Vec3 Vec3::Spherical(float phi, float cosTheta)
{
    auto sinTheta = SafeSqrt(1.0f - cosTheta * cosTheta);

    return Vec3::Spherical(phi, cosTheta, sinTheta);
}

template <class V>
inline Vec3 Vec3::FromCoords(V const& v) {
    return { static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec3 const& u, Vec3 const& v) {
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

inline bool operator!=(Vec3 const& u, Vec3 const& v) {
    return !(u == v);
}

inline bool operator<(Vec3 const& u, Vec3 const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;
    if (u.y < v.y) return true;
    if (v.y < u.y) return false;

    return u.z < v.z;
}

inline bool operator>(Vec3 const& u, Vec3 const& v) {
    return v < u;
}

inline bool operator<=(Vec3 const& u, Vec3 const& v) {
    return !(v < u);
}

inline bool operator>=(Vec3 const& u, Vec3 const& v) {
    return !(u < v);
}

inline Vec3b CmpEQ(Vec3 const& u, Vec3 const& v) { return { u.x == v.x, u.y == v.y, u.z == v.z }; }
inline Vec3b CmpNE(Vec3 const& u, Vec3 const& v) { return { u.x != v.x, u.y != v.y, u.z != v.z }; }
inline Vec3b CmpLT(Vec3 const& u, Vec3 const& v) { return { u.x <  v.x, u.y <  v.y, u.z <  v.z }; }
inline Vec3b CmpLE(Vec3 const& u, Vec3 const& v) { return { u.x <= v.x, u.y <= v.y, u.z <= v.z }; }
inline Vec3b CmpGT(Vec3 const& u, Vec3 const& v) { return { u.x >  v.x, u.y >  v.y, u.z >  v.z }; }
inline Vec3b CmpGE(Vec3 const& u, Vec3 const& v) { return { u.x >= v.x, u.y >= v.y, u.z >= v.z }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec3 operator+(Vec3 const& u) { return u; }
inline Vec3 operator-(Vec3 const& u) { return { -u.x, -u.y, -u.z }; }

inline Vec3 operator+(Vec3 const& u, Vec3 const& v) { return { u.x + v.x, u.y + v.y, u.z + v.z }; }
inline Vec3 operator-(Vec3 const& u, Vec3 const& v) { return { u.x - v.x, u.y - v.y, u.z - v.z }; }
inline Vec3 operator*(Vec3 const& u, Vec3 const& v) { return { u.x * v.x, u.y * v.y, u.z * v.z }; }
inline Vec3 operator/(Vec3 const& u, Vec3 const& v) { return { u.x / v.x, u.y / v.y, u.z / v.z }; }

inline Vec3& operator+=(Vec3& u, Vec3 const& v) { u = u + v; return u; }
inline Vec3& operator-=(Vec3& u, Vec3 const& v) { u = u - v; return u; }
inline Vec3& operator*=(Vec3& u, Vec3 const& v) { u = u * v; return u; }
inline Vec3& operator/=(Vec3& u, Vec3 const& v) { u = u / v; return u; }

inline Vec3 operator+(float u, Vec3 const& v) { return { u + v.x, u + v.y, u + v.z }; }
inline Vec3 operator-(float u, Vec3 const& v) { return { u - v.x, u - v.y, u - v.z }; }
inline Vec3 operator*(float u, Vec3 const& v) { return { u * v.x, u * v.y, u * v.z }; }
inline Vec3 operator/(float u, Vec3 const& v) { return { u / v.x, u / v.y, u / v.z }; }

inline Vec3 operator+(Vec3 const& u, float v) { return { u.x + v, u.y + v, u.z + v }; }
inline Vec3 operator-(Vec3 const& u, float v) { return { u.x - v, u.y - v, u.z - v }; }
inline Vec3 operator*(Vec3 const& u, float v) { return { u.x * v, u.y * v, u.z * v }; }
inline Vec3 operator/(Vec3 const& u, float v) { return { u.x / v, u.y / v, u.z / v }; }

inline Vec3& operator+=(Vec3& u, float v) { u = u + v; return u; }
inline Vec3& operator-=(Vec3& u, float v) { u = u - v; return u; }
inline Vec3& operator*=(Vec3& u, float v) { u = u * v; return u; }
inline Vec3& operator/=(Vec3& u, float v) { u = u / v; return u; }

//------------------------------------------------------------------------------
// Conditional
//

inline Vec3 Select(Vec3b const& condition, Vec3 const& true_, Vec3 const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y,
             condition.z ? true_.z : false_.z };
}

//==============================================================================
// Vec3i
//==============================================================================

inline Vec3i::Vec3i(Vec2 const& v, int z)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(z)
{
}

inline Vec3i::Vec3i(Vec2i const& v, int z)
    : x(v.x)
    , y(v.y)
    , z(z)
{
}

inline Vec3i::Vec3i(Vec3 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(static_cast<int>(v.z))
{
}

inline Vec3i::Vec3i(Vec4 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(static_cast<int>(v.z))
{
}

inline Vec3i::Vec3i(Vec4i const& v)
    : x(v.x)
    , y(v.y)
    , z(v.z)
{
}

template <class V>
inline Vec3i Vec3i::FromCoords(V const& v) {
    return { static_cast<int>(v.x), static_cast<int>(v.y), static_cast<int>(v.z) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec3i const& u, Vec3i const& v) {
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

inline bool operator!=(Vec3i const& u, Vec3i const& v) {
    return !(u == v);
}

inline bool operator<(Vec3i const& u, Vec3i const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;
    if (u.y < v.y) return true;
    if (v.y < u.y) return false;

    return u.z < v.z;
}

inline bool operator>(Vec3i const& u, Vec3i const& v) {
    return v < u;
}

inline bool operator<=(Vec3i const& u, Vec3i const& v) {
    return !(v < u);
}

inline bool operator>=(Vec3i const& u, Vec3i const& v) {
    return !(u < v);
}

inline Vec3b CmpEQ(Vec3i const& u, Vec3i const& v) { return { u.x == v.x, u.y == v.y, u.z == v.z }; }
inline Vec3b CmpNE(Vec3i const& u, Vec3i const& v) { return { u.x != v.x, u.y != v.y, u.z != v.z }; }
inline Vec3b CmpLT(Vec3i const& u, Vec3i const& v) { return { u.x <  v.x, u.y <  v.y, u.z <  v.z }; }
inline Vec3b CmpLE(Vec3i const& u, Vec3i const& v) { return { u.x <= v.x, u.y <= v.y, u.z <= v.z }; }
inline Vec3b CmpGT(Vec3i const& u, Vec3i const& v) { return { u.x >  v.x, u.y >  v.y, u.z >  v.z }; }
inline Vec3b CmpGE(Vec3i const& u, Vec3i const& v) { return { u.x >= v.x, u.y >= v.y, u.z >= v.z }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec3i operator+(Vec3i const& u) { return u; }
inline Vec3i operator-(Vec3i const& u) { return { -u.x, -u.y, -u.z }; }
inline Vec3i operator!(Vec3i const& u) { return { !u.x, !u.y, !u.z }; }
inline Vec3i operator~(Vec3i const& u) { return { ~u.x, ~u.y, ~u.z }; }

inline Vec3i operator+ (Vec3i const& u, Vec3i const& v) { return { u.x +  v.x, u.y +  v.y, u.z +  v.z }; }
inline Vec3i operator- (Vec3i const& u, Vec3i const& v) { return { u.x -  v.x, u.y -  v.y, u.z -  v.z }; }
inline Vec3i operator* (Vec3i const& u, Vec3i const& v) { return { u.x *  v.x, u.y *  v.y, u.z *  v.z }; }
inline Vec3i operator/ (Vec3i const& u, Vec3i const& v) { return { u.x /  v.x, u.y /  v.y, u.z /  v.z }; }
inline Vec3i operator% (Vec3i const& u, Vec3i const& v) { return { u.x %  v.x, u.y %  v.y, u.z %  v.z }; }
inline Vec3i operator& (Vec3i const& u, Vec3i const& v) { return { u.x &  v.x, u.y &  v.y, u.z &  v.z }; }
inline Vec3i operator| (Vec3i const& u, Vec3i const& v) { return { u.x |  v.x, u.y |  v.y, u.z |  v.z }; }
inline Vec3i operator^ (Vec3i const& u, Vec3i const& v) { return { u.x ^  v.x, u.y ^  v.y, u.z ^  v.z }; }
inline Vec3i operator<<(Vec3i const& u, Vec3i const& v) { return { u.x << v.x, u.y << v.y, u.z << v.z }; }
inline Vec3i operator>>(Vec3i const& u, Vec3i const& v) { return { u.x >> v.x, u.y >> v.y, u.z >> v.z }; }

inline Vec3i& operator+= (Vec3i& u, Vec3i const& v) { u = u +  v; return u; }
inline Vec3i& operator-= (Vec3i& u, Vec3i const& v) { u = u -  v; return u; }
inline Vec3i& operator*= (Vec3i& u, Vec3i const& v) { u = u *  v; return u; }
inline Vec3i& operator/= (Vec3i& u, Vec3i const& v) { u = u /  v; return u; }
inline Vec3i& operator%= (Vec3i& u, Vec3i const& v) { u = u %  v; return u; }
inline Vec3i& operator&= (Vec3i& u, Vec3i const& v) { u = u &  v; return u; }
inline Vec3i& operator|= (Vec3i& u, Vec3i const& v) { u = u |  v; return u; }
inline Vec3i& operator^= (Vec3i& u, Vec3i const& v) { u = u ^  v; return u; }
inline Vec3i& operator<<=(Vec3i& u, Vec3i const& v) { u = u << v; return u; }
inline Vec3i& operator>>=(Vec3i& u, Vec3i const& v) { u = u >> v; return u; }

inline Vec3i operator+ (int u, Vec3i const& v) { return { u +  v.x, u +  v.y, u +  v.z }; }
inline Vec3i operator- (int u, Vec3i const& v) { return { u -  v.x, u -  v.y, u -  v.z }; }
inline Vec3i operator* (int u, Vec3i const& v) { return { u *  v.x, u *  v.y, u *  v.z }; }
inline Vec3i operator/ (int u, Vec3i const& v) { return { u /  v.x, u /  v.y, u /  v.z }; }
inline Vec3i operator% (int u, Vec3i const& v) { return { u %  v.x, u %  v.y, u %  v.z }; }
inline Vec3i operator& (int u, Vec3i const& v) { return { u &  v.x, u &  v.y, u &  v.z }; }
inline Vec3i operator| (int u, Vec3i const& v) { return { u |  v.x, u |  v.y, u |  v.z }; }
inline Vec3i operator^ (int u, Vec3i const& v) { return { u ^  v.x, u ^  v.y, u ^  v.z }; }
inline Vec3i operator<<(int u, Vec3i const& v) { return { u << v.x, u << v.y, u << v.z }; }
inline Vec3i operator>>(int u, Vec3i const& v) { return { u >> v.x, u >> v.y, u >> v.z }; }

inline Vec3i operator+ (Vec3i const& u, int v) { return { u.x +  v, u.y +  v, u.z +  v }; }
inline Vec3i operator- (Vec3i const& u, int v) { return { u.x -  v, u.y -  v, u.z -  v }; }
inline Vec3i operator* (Vec3i const& u, int v) { return { u.x *  v, u.y *  v, u.z *  v }; }
inline Vec3i operator/ (Vec3i const& u, int v) { return { u.x /  v, u.y /  v, u.z /  v }; }
inline Vec3i operator% (Vec3i const& u, int v) { return { u.x %  v, u.y %  v, u.z %  v }; }
inline Vec3i operator& (Vec3i const& u, int v) { return { u.x &  v, u.y &  v, u.z &  v }; }
inline Vec3i operator| (Vec3i const& u, int v) { return { u.x |  v, u.y |  v, u.z |  v }; }
inline Vec3i operator^ (Vec3i const& u, int v) { return { u.x ^  v, u.y ^  v, u.z ^  v }; }
inline Vec3i operator<<(Vec3i const& u, int v) { return { u.x << v, u.y << v, u.z << v }; }
inline Vec3i operator>>(Vec3i const& u, int v) { return { u.x >> v, u.y >> v, u.z >> v }; }

inline Vec3i& operator+= (Vec3i& u, int v) { u = u +  v; return u; }
inline Vec3i& operator-= (Vec3i& u, int v) { u = u -  v; return u; }
inline Vec3i& operator*= (Vec3i& u, int v) { u = u *  v; return u; }
inline Vec3i& operator/= (Vec3i& u, int v) { u = u /  v; return u; }
inline Vec3i& operator%= (Vec3i& u, int v) { u = u %  v; return u; }
inline Vec3i& operator&= (Vec3i& u, int v) { u = u &  v; return u; }
inline Vec3i& operator|= (Vec3i& u, int v) { u = u |  v; return u; }
inline Vec3i& operator^= (Vec3i& u, int v) { u = u ^  v; return u; }
inline Vec3i& operator<<=(Vec3i& u, int v) { u = u << v; return u; }
inline Vec3i& operator>>=(Vec3i& u, int v) { u = u >> v; return u; }

inline Vec3i Select(Vec3b const& condition, Vec3i const& true_, Vec3i const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y,
             condition.z ? true_.z : false_.z };
}

//==============================================================================
// Vec4b
//==============================================================================

inline bool All (Vec4b const& v) { return v.x && v.y && v.z && v.w; }
inline bool Any (Vec4b const& v) { return v.x || v.y || v.z || v.w; }
inline bool None(Vec4b const& v) { return !Any(v); }

inline Vec4b operator&(Vec4b const& u, Vec4b const& v) { return { u.x && v.x, u.y && v.y, u.z && v.z, u.w && v.w }; }
inline Vec4b operator|(Vec4b const& u, Vec4b const& v) { return { u.x || v.x, u.y || v.y, u.z || v.z, u.w || v.w }; }

//==============================================================================
// Vec4
//==============================================================================

inline Vec4::Vec4(Vec2 const& v, float z, float w)
    : x(v.x)
    , y(v.y)
    , z(z)
    , w(w)
{
}

inline Vec4::Vec4(Vec2i const& v, float z, float w)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(z)
    , w(w)
{
}

inline Vec4::Vec4(Vec3 const& v, float w)
    : x(v.x)
    , y(v.y)
    , z(v.z)
    , w(w)
{
}

inline Vec4::Vec4(Vec3i const& v, float w)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(static_cast<float>(v.z))
    , w(w)
{
}

inline Vec4::Vec4(Vec4i const& v)
    : x(static_cast<float>(v.x))
    , y(static_cast<float>(v.y))
    , z(static_cast<float>(v.z))
    , w(static_cast<float>(v.w))
{
}

inline Vec4::Vec4(Point3 const& p)
    : x(p.x)
    , y(p.y)
    , z(p.z)
    , w(1)
{
}

template <class V>
inline Vec4 Vec4::FromCoords(V const& v) {
    return { static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z), static_cast<float>(v.w) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec4 const& u, Vec4 const& v) {
    return u.x == v.x && u.y == v.y && u.z == v.z && u.w == v.w;
}

inline bool operator!=(Vec4 const& u, Vec4 const& v) {
    return !(u == v);
}

inline bool operator<(Vec4 const& u, Vec4 const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;
    if (u.y < v.y) return true;
    if (v.y < u.y) return false;
    if (u.z < v.z) return true;
    if (v.z < u.z) return false;

    return u.w < v.w;
}

inline bool operator>(Vec4 const& u, Vec4 const& v) {
    return v < u;
}

inline bool operator<=(Vec4 const& u, Vec4 const& v) {
    return !(v < u);
}

inline bool operator>=(Vec4 const& u, Vec4 const& v) {
    return !(u < v);
}

inline Vec4b CmpEQ(Vec4 const& u, Vec4 const& v) { return { u.x == v.x, u.y == v.y, u.z == v.z, u.w == v.w }; }
inline Vec4b CmpNE(Vec4 const& u, Vec4 const& v) { return { u.x != v.x, u.y != v.y, u.z != v.z, u.w != v.w }; }
inline Vec4b CmpLT(Vec4 const& u, Vec4 const& v) { return { u.x <  v.x, u.y <  v.y, u.z <  v.z, u.w <  v.w }; }
inline Vec4b CmpLE(Vec4 const& u, Vec4 const& v) { return { u.x <= v.x, u.y <= v.y, u.z <= v.z, u.w <= v.w }; }
inline Vec4b CmpGT(Vec4 const& u, Vec4 const& v) { return { u.x >  v.x, u.y >  v.y, u.z >  v.z, u.w >  v.w }; }
inline Vec4b CmpGE(Vec4 const& u, Vec4 const& v) { return { u.x >= v.x, u.y >= v.y, u.z >= v.z, u.w >= v.w }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec4 operator+(Vec4 const& u) { return u; }
inline Vec4 operator-(Vec4 const& u) { return { -u.x, -u.y, -u.z, -u.w }; }

inline Vec4 operator+(Vec4 const& u, Vec4 const& v) { return { u.x + v.x, u.y + v.y, u.z + v.z, u.w + v.w }; }
inline Vec4 operator-(Vec4 const& u, Vec4 const& v) { return { u.x - v.x, u.y - v.y, u.z - v.z, u.w - v.w }; }
inline Vec4 operator*(Vec4 const& u, Vec4 const& v) { return { u.x * v.x, u.y * v.y, u.z * v.z, u.w * v.w }; }
inline Vec4 operator/(Vec4 const& u, Vec4 const& v) { return { u.x / v.x, u.y / v.y, u.z / v.z, u.w / v.w }; }

inline Vec4& operator+=(Vec4& u, Vec4 const& v) { u = u + v; return u; }
inline Vec4& operator-=(Vec4& u, Vec4 const& v) { u = u - v; return u; }
inline Vec4& operator*=(Vec4& u, Vec4 const& v) { u = u * v; return u; }
inline Vec4& operator/=(Vec4& u, Vec4 const& v) { u = u / v; return u; }

inline Vec4 operator+(float u, Vec4 const& v) { return { u + v.x, u + v.y, u + v.z, u + v.w }; }
inline Vec4 operator-(float u, Vec4 const& v) { return { u - v.x, u - v.y, u - v.z, u - v.w }; }
inline Vec4 operator*(float u, Vec4 const& v) { return { u * v.x, u * v.y, u * v.z, u * v.w }; }
inline Vec4 operator/(float u, Vec4 const& v) { return { u / v.x, u / v.y, u / v.z, u / v.w }; }

inline Vec4 operator+(Vec4 const& u, float v) { return { u.x + v, u.y + v, u.z + v, u.w + v }; }
inline Vec4 operator-(Vec4 const& u, float v) { return { u.x - v, u.y - v, u.z - v, u.w - v }; }
inline Vec4 operator*(Vec4 const& u, float v) { return { u.x * v, u.y * v, u.z * v, u.w * v }; }
inline Vec4 operator/(Vec4 const& u, float v) { return { u.x / v, u.y / v, u.z / v, u.w / v }; }

inline Vec4& operator+=(Vec4& u, float v) { u = u + v; return u; }
inline Vec4& operator-=(Vec4& u, float v) { u = u - v; return u; }
inline Vec4& operator*=(Vec4& u, float v) { u = u * v; return u; }
inline Vec4& operator/=(Vec4& u, float v) { u = u / v; return u; }

//------------------------------------------------------------------------------
// Conditional
//

inline Vec4 Select(Vec4b const& condition, Vec4 const& true_, Vec4 const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y,
             condition.z ? true_.z : false_.z,
             condition.w ? true_.w : false_.w };
}

//==============================================================================
// Vec4i
//==============================================================================

inline Vec4i::Vec4i(Vec2 const& v, int z, int w)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(z)
    , w(w)
{
}

inline Vec4i::Vec4i(Vec2i const& v, int z, int w)
    : x(v.x)
    , y(v.y)
    , z(z)
    , w(w)
{
}

inline Vec4i::Vec4i(Vec3 const& v, int w)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(static_cast<int>(v.z))
    , w(w)
{
}

inline Vec4i::Vec4i(Vec3i const& v, int w)
    : x(v.x)
    , y(v.y)
    , z(v.z)
    , w(w)
{
}

inline Vec4i::Vec4i(Vec4 const& v)
    : x(static_cast<int>(v.x))
    , y(static_cast<int>(v.y))
    , z(static_cast<int>(v.z))
    , w(static_cast<int>(v.w))
{
}

template <class V>
inline Vec4i Vec4i::FromCoords(V const& v) {
    return { static_cast<int>(v.x), static_cast<int>(v.y), static_cast<int>(v.z), static_cast<int>(v.w) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Vec4i const& u, Vec4i const& v) {
    return u.x == v.x && u.y == v.y && u.z == v.z && u.w == v.w;
}

inline bool operator!=(Vec4i const& u, Vec4i const& v) {
    return !(u == v);
}

inline bool operator<(Vec4i const& u, Vec4i const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;
    if (u.y < v.y) return true;
    if (v.y < u.y) return false;
    if (u.z < v.z) return true;
    if (v.z < u.z) return false;

    return u.w < v.w;
}

inline bool operator>(Vec4i const& u, Vec4i const& v) {
    return v < u;
}

inline bool operator<=(Vec4i const& u, Vec4i const& v) {
    return !(v < u);
}

inline bool operator>=(Vec4i const& u, Vec4i const& v) {
    return !(u < v);
}

inline Vec4b CmpEQ(Vec4i const& u, Vec4i const& v) { return { u.x == v.x, u.y == v.y, u.z == v.z, u.w == v.w }; }
inline Vec4b CmpNE(Vec4i const& u, Vec4i const& v) { return { u.x != v.x, u.y != v.y, u.z != v.z, u.w != v.w }; }
inline Vec4b CmpLT(Vec4i const& u, Vec4i const& v) { return { u.x <  v.x, u.y <  v.y, u.z <  v.z, u.w <  v.w }; }
inline Vec4b CmpLE(Vec4i const& u, Vec4i const& v) { return { u.x <= v.x, u.y <= v.y, u.z <= v.z, u.w <= v.w }; }
inline Vec4b CmpGT(Vec4i const& u, Vec4i const& v) { return { u.x >  v.x, u.y >  v.y, u.z >  v.z, u.w >  v.w }; }
inline Vec4b CmpGE(Vec4i const& u, Vec4i const& v) { return { u.x >= v.x, u.y >= v.y, u.z >= v.z, u.w >= v.w }; }

//------------------------------------------------------------------------------
// Arithmetic
//

inline Vec4i operator+(Vec4i const& u) { return u; }
inline Vec4i operator-(Vec4i const& u) { return { -u.x, -u.y, -u.z, -u.w }; }
inline Vec4i operator!(Vec4i const& u) { return { !u.x, !u.y, !u.z, !u.w }; }
inline Vec4i operator~(Vec4i const& u) { return { ~u.x, ~u.y, ~u.z, ~u.w }; }

inline Vec4i operator+ (Vec4i const& u, Vec4i const& v) { return { u.x +  v.x, u.y +  v.y, u.z +  v.z, u.w +  v.w }; }
inline Vec4i operator- (Vec4i const& u, Vec4i const& v) { return { u.x -  v.x, u.y -  v.y, u.z -  v.z, u.w -  v.w }; }
inline Vec4i operator* (Vec4i const& u, Vec4i const& v) { return { u.x *  v.x, u.y *  v.y, u.z *  v.z, u.w *  v.w }; }
inline Vec4i operator/ (Vec4i const& u, Vec4i const& v) { return { u.x /  v.x, u.y /  v.y, u.z /  v.z, u.w /  v.w }; }
inline Vec4i operator% (Vec4i const& u, Vec4i const& v) { return { u.x %  v.x, u.y %  v.y, u.z %  v.z, u.w %  v.w }; }
inline Vec4i operator& (Vec4i const& u, Vec4i const& v) { return { u.x &  v.x, u.y &  v.y, u.z &  v.z, u.w &  v.w }; }
inline Vec4i operator| (Vec4i const& u, Vec4i const& v) { return { u.x |  v.x, u.y |  v.y, u.z |  v.z, u.w |  v.w }; }
inline Vec4i operator^ (Vec4i const& u, Vec4i const& v) { return { u.x ^  v.x, u.y ^  v.y, u.z ^  v.z, u.w ^  v.w }; }
inline Vec4i operator<<(Vec4i const& u, Vec4i const& v) { return { u.x << v.x, u.y << v.y, u.z << v.z, u.w << v.w }; }
inline Vec4i operator>>(Vec4i const& u, Vec4i const& v) { return { u.x >> v.x, u.y >> v.y, u.z >> v.z, u.w >> v.w }; }

inline Vec4i& operator+= (Vec4i& u, Vec4i const& v) { u = u +  v; return u; }
inline Vec4i& operator-= (Vec4i& u, Vec4i const& v) { u = u -  v; return u; }
inline Vec4i& operator*= (Vec4i& u, Vec4i const& v) { u = u *  v; return u; }
inline Vec4i& operator/= (Vec4i& u, Vec4i const& v) { u = u /  v; return u; }
inline Vec4i& operator%= (Vec4i& u, Vec4i const& v) { u = u %  v; return u; }
inline Vec4i& operator&= (Vec4i& u, Vec4i const& v) { u = u &  v; return u; }
inline Vec4i& operator|= (Vec4i& u, Vec4i const& v) { u = u |  v; return u; }
inline Vec4i& operator^= (Vec4i& u, Vec4i const& v) { u = u ^  v; return u; }
inline Vec4i& operator<<=(Vec4i& u, Vec4i const& v) { u = u << v; return u; }
inline Vec4i& operator>>=(Vec4i& u, Vec4i const& v) { u = u >> v; return u; }

inline Vec4i operator+ (int u, Vec4i const& v) { return { u +  v.x, u +  v.y, u +  v.z, u +  v.w }; }
inline Vec4i operator- (int u, Vec4i const& v) { return { u -  v.x, u -  v.y, u -  v.z, u -  v.w }; }
inline Vec4i operator* (int u, Vec4i const& v) { return { u *  v.x, u *  v.y, u *  v.z, u *  v.w }; }
inline Vec4i operator/ (int u, Vec4i const& v) { return { u /  v.x, u /  v.y, u /  v.z, u /  v.w }; }
inline Vec4i operator% (int u, Vec4i const& v) { return { u %  v.x, u %  v.y, u %  v.z, u %  v.w }; }
inline Vec4i operator& (int u, Vec4i const& v) { return { u &  v.x, u &  v.y, u &  v.z, u &  v.w }; }
inline Vec4i operator| (int u, Vec4i const& v) { return { u |  v.x, u |  v.y, u |  v.z, u |  v.w }; }
inline Vec4i operator^ (int u, Vec4i const& v) { return { u ^  v.x, u ^  v.y, u ^  v.z, u ^  v.w }; }
inline Vec4i operator<<(int u, Vec4i const& v) { return { u << v.x, u << v.y, u << v.z, u << v.w }; }
inline Vec4i operator>>(int u, Vec4i const& v) { return { u >> v.x, u >> v.y, u >> v.z, u >> v.w }; }

inline Vec4i operator+ (Vec4i const& u, int v) { return { u.x +  v, u.y +  v, u.z +  v, u.w +  v }; }
inline Vec4i operator- (Vec4i const& u, int v) { return { u.x -  v, u.y -  v, u.z -  v, u.w -  v }; }
inline Vec4i operator* (Vec4i const& u, int v) { return { u.x *  v, u.y *  v, u.z *  v, u.w *  v }; }
inline Vec4i operator/ (Vec4i const& u, int v) { return { u.x /  v, u.y /  v, u.z /  v, u.w /  v }; }
inline Vec4i operator% (Vec4i const& u, int v) { return { u.x %  v, u.y %  v, u.z %  v, u.w %  v }; }
inline Vec4i operator& (Vec4i const& u, int v) { return { u.x &  v, u.y &  v, u.z &  v, u.w &  v }; }
inline Vec4i operator| (Vec4i const& u, int v) { return { u.x |  v, u.y |  v, u.z |  v, u.w |  v }; }
inline Vec4i operator^ (Vec4i const& u, int v) { return { u.x ^  v, u.y ^  v, u.z ^  v, u.w ^  v }; }
inline Vec4i operator<<(Vec4i const& u, int v) { return { u.x << v, u.y << v, u.z << v, u.w << v }; }
inline Vec4i operator>>(Vec4i const& u, int v) { return { u.x >> v, u.y >> v, u.z >> v, u.w >> v }; }

inline Vec4i& operator+= (Vec4i& u, int v) { u = u +  v; return u; }
inline Vec4i& operator-= (Vec4i& u, int v) { u = u -  v; return u; }
inline Vec4i& operator*= (Vec4i& u, int v) { u = u *  v; return u; }
inline Vec4i& operator/= (Vec4i& u, int v) { u = u /  v; return u; }
inline Vec4i& operator%= (Vec4i& u, int v) { u = u %  v; return u; }
inline Vec4i& operator&= (Vec4i& u, int v) { u = u &  v; return u; }
inline Vec4i& operator|= (Vec4i& u, int v) { u = u |  v; return u; }
inline Vec4i& operator^= (Vec4i& u, int v) { u = u ^  v; return u; }
inline Vec4i& operator<<=(Vec4i& u, int v) { u = u << v; return u; }
inline Vec4i& operator>>=(Vec4i& u, int v) { u = u >> v; return u; }

//------------------------------------------------------------------------------
// Conditional
//

inline Vec4i Select(Vec4b const& condition, Vec4i const& true_, Vec4i const& false_)
{
    return { condition.x ? true_.x : false_.x,
             condition.y ? true_.y : false_.y,
             condition.z ? true_.z : false_.z,
             condition.w ? true_.w : false_.w };
}

//==============================================================================
// Point3
//==============================================================================

inline Point3::Point3(Vec3 const& v)
    : x(v.x)
    , y(v.y)
    , z(v.z)
{
}

inline Point3::Point3(Vec4 const& v)
    : x(v.x / v.w)
    , y(v.y / v.w)
    , z(v.z / v.w)
{
}

template <class V>
inline Point3 Point3::FromCoords(V const& v) {
    return { static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z) };
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Point3 const& u, Point3 const& v) {
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

inline bool operator!=(Point3 const& u, Point3 const& v) {
    return !(u == v);
}

inline bool operator<(Point3 const& u, Point3 const& v)
{
    if (u.x < v.x) return true;
    if (v.x < u.x) return false;
    if (u.y < v.y) return true;
    if (v.y < u.y) return false;

    return u.z < v.z;
}

inline bool operator>(Point3 const& u, Point3 const& v) {
    return v < u;
}

inline bool operator<=(Point3 const& u, Point3 const& v) {
    return !(v < u);
}

inline bool operator>=(Point3 const& u, Point3 const& v) {
    return !(u < v);
}

//------------------------------------------------------------------------------
// Arithmetic
//

inline Point3 operator+(Point3 const& u, Vec3   const& v) { return { u.x + v.x, u.y + v.y, u.z + v.z }; }
inline Point3 operator+(Vec3   const& u, Point3 const& v) { return { u.x + v.x, u.y + v.y, u.z + v.z }; }
inline Vec3   operator-(Point3 const& u, Point3 const& v) { return { u.x - v.x, u.y - v.y, u.z - v.z }; }
inline Point3 operator-(Point3 const& u, Vec3   const& v) { return { u.x - v.x, u.y - v.y, u.z - v.z }; }

inline Point3& operator +=(Point3& u, Vec3 const& v) { u = u + v; return u; }
inline Point3& operator -=(Point3& u, Vec3 const& v) { u = u - v; return u; }

//==============================================================================
// Geometric functions
//==============================================================================

//------------------------------------------------------------------------------
// Dot
//

inline float Dot(Vec2  const& u, Vec2  const& v) { return u.x * v.x + u.y * v.y; }
inline int   Dot(Vec2i const& u, Vec2i const& v) { return u.x * v.x + u.y * v.y; }
inline float Dot(Vec3  const& u, Vec3  const& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
inline int   Dot(Vec3i const& u, Vec3i const& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
inline float Dot(Vec4  const& u, Vec4  const& v) { return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w; }
inline int   Dot(Vec4i const& u, Vec4i const& v) { return u.x * v.x + u.y * v.y + u.z * v.z + u.w * v.w; }

inline float Dot(Point3 const& u, Vec3   const& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
inline float Dot(Vec3   const& u, Point3 const& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }

//------------------------------------------------------------------------------
// SaturateDot
//

inline float SaturateDot(Vec2 const& u, Vec2 const& v) { return Saturate(Dot(u, v)); }
inline float SaturateDot(Vec3 const& u, Vec3 const& v) { return Saturate(Dot(u, v)); }
inline float SaturateDot(Vec4 const& u, Vec4 const& v) { return Saturate(Dot(u, v)); }

//------------------------------------------------------------------------------
// LengthSquared
//

inline float LengthSquared(Vec2  const& u) { return Dot(u, u); }
inline int   LengthSquared(Vec2i const& u) { return Dot(u, u); }
inline float LengthSquared(Vec3  const& u) { return Dot(u, u); }
inline int   LengthSquared(Vec3i const& u) { return Dot(u, u); }
inline float LengthSquared(Vec4  const& u) { return Dot(u, u); }
inline int   LengthSquared(Vec4i const& u) { return Dot(u, u); }

//------------------------------------------------------------------------------
// Length
//

inline float LengthSlow(Vec2 const& u) { return Hypot(u.x, u.y); }
inline float LengthSlow(Vec3 const& u) { return Hypot(u.x, u.y, u.z); }
inline float LengthSlow(Vec4 const& u) { return Hypot(u.x, u.y, u.z, u.w); }

inline float LengthFast(Vec2 const& u) { return Sqrt(LengthSquared(u)); }
inline float LengthFast(Vec3 const& u) { return Sqrt(LengthSquared(u)); }
inline float LengthFast(Vec4 const& u) { return Sqrt(LengthSquared(u)); }

#if 1
inline float Length(Vec2 const& u) { return LengthFast(u); }
inline float Length(Vec3 const& u) { return LengthFast(u); }
inline float Length(Vec4 const& u) { return LengthFast(u); }
#else
inline float Length(Vec2 const& u) { return LengthSlow(u); }
inline float Length(Vec3 const& u) { return LengthSlow(u); }
inline float Length(Vec4 const& u) { return LengthSlow(u); }
#endif

//------------------------------------------------------------------------------
// Normalize
//

inline Vec2 Normalize(Vec2 const& u) { return u / Length(u); }
inline Vec3 Normalize(Vec3 const& u) { return u / Length(u); }
inline Vec4 Normalize(Vec4 const& u) { return u / Length(u); }

//------------------------------------------------------------------------------
// Cross
//
// Returns the cross product of U and V.
//

//inline float Cross(Vec2 const& u, Vec2 const& v) {
//    return u.x * v.y - u.y * v.x;
//}

inline Vec3 Cross(Vec3 const& u, Vec3 const& v)
{
    const float x = u.y * v.z - u.z * v.y;
    const float y = u.z * v.x - u.x * v.z;
    const float z = u.x * v.y - u.y * v.x;

    return {x, y, z};
}

//------------------------------------------------------------------------------
// Angle
//
// Note:
// u and v must be normalized!
//

#if 0

inline float Angle(Vec2 const& u, Vec2 const& v) { return 2 * Atan2(Length(u - v), Length(u + v)); }
inline float Angle(Vec3 const& u, Vec3 const& v) { return 2 * Atan2(Length(u - v), Length(u + v)); }
inline float Angle(Vec4 const& u, Vec4 const& v) { return 2 * Atan2(Length(u - v), Length(u + v)); }

#else

inline float Angle(Vec2 const& u, Vec2 const& v)
{
    const float sin_phi = u.x * v.y - u.y * v.x;
    const float cos_phi = Dot(u,v);

    return Atan2(sin_phi, cos_phi);
}

inline float Angle(Vec3 const& u, Vec3 const& v)
{
    const float sin_phi = Length(Cross(u,v));
    const float cos_phi = Dot(u,v);

    return Atan2(sin_phi, cos_phi);
}

inline float Angle(Vec4 const& u, Vec4 const& v) { return 2 * Atan2(Length(u - v), Length(u + v)); }

#endif

//------------------------------------------------------------------------------
// FaceForward
//
// Return a vector pointing in the same direction as another
//

inline Vec2 FaceForward(Vec2 const& Nref, Vec2 const& I, Vec2 const& N) { return Dot(Nref, I) < 0 ? N : -N; }
inline Vec3 FaceForward(Vec3 const& Nref, Vec3 const& I, Vec3 const& N) { return Dot(Nref, I) < 0 ? N : -N; }

inline Vec2 FaceForward(Vec2 const& I, Vec2 const& N) { return FaceForward(N, I, N); }
inline Vec3 FaceForward(Vec3 const& I, Vec3 const& N) { return FaceForward(N, I, N); }

//------------------------------------------------------------------------------
// Reflect
//
// Calculate the reflection direction for an incident vector
//

inline Vec2 Reflect(Vec2 const& I, Vec2 const& N) { return I - (2 * Dot(I, N)) * N; }
inline Vec3 Reflect(Vec3 const& I, Vec3 const& N) { return I - (2 * Dot(I, N)) * N; }

//------------------------------------------------------------------------------
// Perp
//
// Returns a vector orthogonal to the given vector n.
//

inline Vec2 Perp(Vec2 const& n) {
    return { -n.y, n.x };
}

inline Vec3 Perp(Vec3 const& n)
{
#if 1
    // From:
    // http://blog.selfshadow.com/2011/10/17/perp-vectors/
    //
    // v = n x (0, A, 1) where A = -Sign(yz)

    const float A = -SignNotZero(n.y * n.z);
    return { n.y - A * n.z, -n.x, n.x * A };
#else
    // From:
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts

    return Abs(n.x) > Abs(n.y) ? Vec3(-n.y, n.x, 0.0f) : Vec3(0.0f, -n.z, n.y);
#endif
}

//==============================================================================
// Misc. functions
//==============================================================================

//------------------------------------------------------------------------------
// MinIndex / MaxIndex
//

inline int MinIndex(Vec2   const& u) { return MinIndex(u.x, u.y); }
inline int MinIndex(Vec2i  const& u) { return MinIndex(u.x, u.y); }
inline int MinIndex(Vec3   const& u) { return MinIndex(u.x, u.y, u.z); }
inline int MinIndex(Vec3i  const& u) { return MinIndex(u.x, u.y, u.z); }
inline int MinIndex(Vec4   const& u) { return MinIndex(u.x, u.y, u.z, u.w); }
inline int MinIndex(Vec4i  const& u) { return MinIndex(u.x, u.y, u.z, u.w); }
inline int MinIndex(Point3 const& u) { return MinIndex(u.x, u.y, u.z); }

inline int MaxIndex(Vec2   const& u) { return MaxIndex(u.x, u.y); }
inline int MaxIndex(Vec2i  const& u) { return MaxIndex(u.x, u.y); }
inline int MaxIndex(Vec3   const& u) { return MaxIndex(u.x, u.y, u.z); }
inline int MaxIndex(Vec3i  const& u) { return MaxIndex(u.x, u.y, u.z); }
inline int MaxIndex(Vec4   const& u) { return MaxIndex(u.x, u.y, u.z, u.w); }
inline int MaxIndex(Vec4i  const& u) { return MaxIndex(u.x, u.y, u.z, u.w); }
inline int MaxIndex(Point3 const& u) { return MaxIndex(u.x, u.y, u.z); }

//------------------------------------------------------------------------------
// Reductions
//

inline float MinElement(Vec2   const& u) { return Min(u.x, u.y); }
inline int   MinElement(Vec2i  const& u) { return Min(u.x, u.y); }
inline float MinElement(Vec3   const& u) { return Min(u.x, u.y, u.z); }
inline int   MinElement(Vec3i  const& u) { return Min(u.x, u.y, u.z); }
inline float MinElement(Vec4   const& u) { return Min(u.x, u.y, u.z, u.w); }
inline int   MinElement(Vec4i  const& u) { return Min(u.x, u.y, u.z, u.w); }
inline float MinElement(Point3 const& u) { return Min(u.x, u.y, u.z); }

inline float MaxElement(Vec2   const& u) { return Max(u.x, u.y); }
inline int   MaxElement(Vec2i  const& u) { return Max(u.x, u.y); }
inline float MaxElement(Vec3   const& u) { return Max(u.x, u.y, u.z); }
inline int   MaxElement(Vec3i  const& u) { return Max(u.x, u.y, u.z); }
inline float MaxElement(Vec4   const& u) { return Max(u.x, u.y, u.z, u.w); }
inline int   MaxElement(Vec4i  const& u) { return Max(u.x, u.y, u.z, u.w); }
inline float MaxElement(Point3 const& u) { return Max(u.x, u.y, u.z); }

inline float ReduceAdd(Vec2  const& u) { return u.x + u.y; }
inline int   ReduceAdd(Vec2i const& u) { return u.x + u.y; }
inline float ReduceAdd(Vec3  const& u) { return u.x + u.y + u.z; }
inline int   ReduceAdd(Vec3i const& u) { return u.x + u.y + u.z; }
inline float ReduceAdd(Vec4  const& u) { return u.x + u.y + u.z + u.w; }
inline int   ReduceAdd(Vec4i const& u) { return u.x + u.y + u.z + u.w; }

inline float ReduceMul(Vec2  const& u) { return u.x * u.y; }
inline int   ReduceMul(Vec2i const& u) { return u.x * u.y; }
inline float ReduceMul(Vec3  const& u) { return u.x * u.y * u.z; }
inline int   ReduceMul(Vec3i const& u) { return u.x * u.y * u.z; }
inline float ReduceMul(Vec4  const& u) { return u.x * u.y * u.z * u.w; }
inline int   ReduceMul(Vec4i const& u) { return u.x * u.y * u.z * u.w; }

//------------------------------------------------------------------------------
// Min / Max
//

inline Vec2  Min(Vec2  const& u, Vec2  const& v) { return Vec2 (Min(u.x, v.x), Min(u.y, v.y)); }
inline Vec2i Min(Vec2i const& u, Vec2i const& v) { return Vec2i(Min(u.x, v.x), Min(u.y, v.y)); }
inline Vec3  Min(Vec3  const& u, Vec3  const& v) { return Vec3 (Min(u.x, v.x), Min(u.y, v.y), Min(u.z, v.z)); }
inline Vec3i Min(Vec3i const& u, Vec3i const& v) { return Vec3i(Min(u.x, v.x), Min(u.y, v.y), Min(u.z, v.z)); }
inline Vec4  Min(Vec4  const& u, Vec4  const& v) { return Vec4 (Min(u.x, v.x), Min(u.y, v.y), Min(u.z, v.z), Min(u.w, v.w)); }
inline Vec4i Min(Vec4i const& u, Vec4i const& v) { return Vec4i(Min(u.x, v.x), Min(u.y, v.y), Min(u.z, v.z), Min(u.w, v.w)); }

inline Vec2  Max(Vec2  const& u, Vec2  const& v) { return Vec2 (Max(u.x, v.x), Max(u.y, v.y)); }
inline Vec2i Max(Vec2i const& u, Vec2i const& v) { return Vec2i(Max(u.x, v.x), Max(u.y, v.y)); }
inline Vec3  Max(Vec3  const& u, Vec3  const& v) { return Vec3 (Max(u.x, v.x), Max(u.y, v.y), Max(u.z, v.z)); }
inline Vec3i Max(Vec3i const& u, Vec3i const& v) { return Vec3i(Max(u.x, v.x), Max(u.y, v.y), Max(u.z, v.z)); }
inline Vec4  Max(Vec4  const& u, Vec4  const& v) { return Vec4 (Max(u.x, v.x), Max(u.y, v.y), Max(u.z, v.z), Max(u.w, v.w)); }
inline Vec4i Max(Vec4i const& u, Vec4i const& v) { return Vec4i(Max(u.x, v.x), Max(u.y, v.y), Max(u.z, v.z), Max(u.w, v.w)); }

//------------------------------------------------------------------------------
// Saturate
//

inline Vec2 Saturate(Vec2 const& u) { return Vec2(Saturate(u.x), Saturate(u.y)); }
inline Vec3 Saturate(Vec3 const& u) { return Vec3(Saturate(u.x), Saturate(u.y), Saturate(u.z)); }
inline Vec4 Saturate(Vec4 const& u) { return Vec4(Saturate(u.x), Saturate(u.y), Saturate(u.z), Saturate(u.w)); }

//------------------------------------------------------------------------------
// Rcp
//

inline Vec2 Rcp(Vec2 const& u) { return Vec2(Rcp(u.x), Rcp(u.y)); }
inline Vec3 Rcp(Vec3 const& u) { return Vec3(Rcp(u.x), Rcp(u.y), Rcp(u.z)); }
inline Vec4 Rcp(Vec4 const& u) { return Vec4(Rcp(u.x), Rcp(u.y), Rcp(u.z), Rcp(u.w)); }

//------------------------------------------------------------------------------
// Sqrt
//

inline Vec2 Sqrt(Vec2 const& u) { return Vec2(Sqrt(u.x), Sqrt(u.y)); }
inline Vec3 Sqrt(Vec3 const& u) { return Vec3(Sqrt(u.x), Sqrt(u.y), Sqrt(u.z)); }
inline Vec4 Sqrt(Vec4 const& u) { return Vec4(Sqrt(u.x), Sqrt(u.y), Sqrt(u.z), Sqrt(u.w)); }

} // namespace math
