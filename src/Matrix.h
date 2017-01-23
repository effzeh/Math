#pragma once

#define LIB_MATH_MATRIX_H 1

#include "Vector.h"

namespace math {

//==============================================================================
// Class definitions
//==============================================================================

struct Euler; // YXZ
struct Mat3;
struct Mat3x4;
struct Mat4;
struct Quat;

struct Mat3
{
    // [ L 0 ]
    // [ 0 1 ]

    Vec3 col0; // x
    Vec3 col1; // y
    Vec3 col2; // z (normal)

    Mat3() = default;

    Mat3(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2);

    Mat3(float m00, float m10, float m20,
         float m01, float m11, float m21,
         float m02, float m12, float m22);

    Mat3(float const* p);

    // Construct from diagonal
    Mat3(float m00, float m11, float m22);

    // Construct from affine transformation (i.e. extract the rotational part)
    explicit Mat3(Mat3x4 const& A);

    // Construct from a Mat4 (i.e. extract the rotational part)
    explicit Mat3(Mat4 const& A);

    float* data();

    float const* data() const;

    Vec3& operator()(int col);

    Vec3 const& operator()(int col) const;

    float& operator()(int row, int col);

    float const& operator()(int row, int col) const;

    // Construct from rows
    static Mat3 FromRows(Vec3 const& r0, Vec3 const& r1, Vec3 const& r2);
    // Returns the identity transformation
    static Mat3 Identity();
    // Returns the zero matrix
    static Mat3 Zero();
    // Returns a scaling transformation
    static Mat3 Scaling(Vec3 const& s);
    // Construct a shear transform
    static Mat3 Shear(float xy, float yx, float xz, float zx, float yz, float zy);
    // Construct a rotation
    static Mat3 Rotation(Vec3 const& axis, float angle);
    // Construct a rotation
    static Mat3 Rotation(Vec3 const& from, Vec3 const& to);
    // Constructs a rotation matrix M from a quaternion q such that M v = q v q'.
    static Mat3 Rotation(Quat const& q);
    // Constructs a rotation matrix from Euler angles
    static Mat3 Rotation(Euler const& e);
    // Constructs a linear map R from a normal N such that R*(0,0,1) = N.
    static Mat3 Frame(Vec3 const& n);

    // Convert from world coordinates to local coordinates
    Vec3 ToLocal(Vec3 const& v) const;

    // Convert from local coordinates to world coordinates
    Vec3 ToWorld(Vec3 const& v) const;

    //
    // The following methods assume that the vector v is given in local coordinates.
    //
    // THETA is the angle between the vector v and the z-axis (i.e. the normal).
    // PHI is the angle between the projection of v into the xy-plane and the x-axis.
    //
    static float CosTheta2(Vec3 const& v);
    static float CosTheta(Vec3 const& v);
    static float SinTheta2(Vec3 const& v);
    static float SinTheta(Vec3 const& v);
    static float TanTheta2(Vec3 const& v);
    static float TanTheta(Vec3 const& v);
    static float SinPhi2(Vec3 const& v);
    static float SinPhi(Vec3 const& v);
    static float CosPhi2(Vec3 const& v);
    static float CosPhi(Vec3 const& v);
};

struct Mat3x4
{
    // [ L t ]
    // [ 0 1 ]

    Vec3 col0; // linear transform L
    Vec3 col1;
    Vec3 col2;
    Vec3 col3; // translation t

    Mat3x4() = default;

    Mat3x4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2, Vec3 const& c3);

    Mat3x4(float m00, float m10, float m20,
           float m01, float m11, float m21,
           float m02, float m12, float m22,
           float m03, float m13, float m23);

    explicit Mat3x4(Mat3 const& L);

    explicit Mat3x4(Mat3 const& L, Vec3 const& T);

    explicit Mat3x4(Mat3 const& L, Point3 const& T);

    explicit Mat3x4(Mat4 const& A);

    // Returns the linear part of the affine transformation
    Mat3& linear() { return *reinterpret_cast<Mat3*>(this); }

    // Returns the linear part of the affine transformation
    Mat3 const& linear() const { return *reinterpret_cast<Mat3 const*>(this); }

    float* data();

    float const* data() const;

    Vec3& operator()(int col);

    Vec3 const& operator()(int col) const;

    float& operator()(int row, int col);

    float const& operator()(int row, int col) const;

    // Returns the identity transformation
    static Mat3x4 Identity();
    // Returns the zero matrix
    static Mat3x4 Zero();
    // Construct a translation
    static Mat3x4 Translation(Vec3 const& t);
    // Construct a scaling transformation
    static Mat3x4 Scaling(Vec3 const& s);
    // Construct a scaling relative to the given point
    static Mat3x4 Scaling(Vec3 const& s, Vec3 const& origin);
    // Construct a shear matrix.
    static Mat3x4 Shear(float xy, float yx, float xz, float zx, float yz, float zy);
    // Construct a rotation
    static Mat3x4 Rotation(Vec3 const& axis, float angle);
    // Construct a rotation
    static Mat3x4 Rotation(Vec3 const& from, Vec3 const& to);
    // Constructs a rotation matrix M from a quaternion q such that M v = q v q'.
    static Mat3x4 Rotation(Quat const& q);
    // Construct a rotation matrix from Euler angles.
    static Mat3x4 Rotation(Euler const& e);
    // Construct an OpenGL-compatible look-at matrix
    static Mat3x4 LookAt(Point3 const& position, Point3 const& target, Vec3 const& up);
};

struct Mat4
{
    // [L t]
    // [u w]

    Vec4 col0;
    Vec4 col1;
    Vec4 col2;
    Vec4 col3;

    Mat4() = default;

    Mat4(Vec4 const& c0, Vec4 const& c1, Vec4 const& c2, Vec4 const& c3);

    Mat4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2);

    Mat4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2, Vec3 const& c3);

    Mat4(float m00, float m10, float m20, float m30,
         float m01, float m11, float m21, float m31,
         float m02, float m12, float m22, float m32,
         float m03, float m13, float m23, float m33);

    explicit Mat4(float const* p);

    // Construct from diagonal
    Mat4(float m00, float m11, float m22, float m33);

    explicit Mat4(Mat3 const& L);

    explicit Mat4(Mat3 const& L, Vec3 const& T);

    explicit Mat4(Mat3 const& L, Point3 const& T);

    explicit Mat4(Mat3x4 const& A);

    float* data();

    float const* data() const;

    Vec4& operator()(int n);

    Vec4 const& operator()(int n) const;

    float& operator()(int row, int col);

    float const& operator()(int row, int col) const;

    // Construct from rows
    static Mat4 FromRows(Vec4 const& r0, Vec4 const& r1, Vec4 const& r2, Vec4 const& r3);
    // Returns the identity matrix
    static Mat4 Identity();
    // Returns the zero matrix
    static Mat4 Zero();
    // Constructs a translation matrix
    static Mat4 Translation(Vec3 const& t);
    // Constructs a scaling matrix
    static Mat4 Scaling(Vec3 const& s);
    // Construct a scaling relative to the given point
    static Mat4 Scaling(Vec3 const& s, Vec3 const& origin);
    // Construct a shear matrix.
    static Mat4 Shear(float xy, float yx, float xz, float zx, float yz, float zy);
    // Constructs a rotation matrix from a rotation tvec and an angle
    static Mat4 Rotation(Vec3 const& axis, float angle);
    // Constructs a rotation matrix, which rotates <from> into <to>
    static Mat4 Rotation(Vec3 const& from, Vec3 const& to);
    // Constructs a rotation matrix M from a quaternion q such that M v = q v q'.
    static Mat4 Rotation(Quat const& q);
    // Construct a rotation matrix from Euler angles
    static Mat4 Rotation(Euler const& q);
    // Constructs a look-at matrix (compatible with OpenGL)
    static Mat4 LookAt(Point3 const& position, Point3 const& target, Vec3 const& up);
    // Constructs an orthogonal projection matrix (compatible with OpenGL)
    static Mat4 Ortho(float left, float right, float bottom, float top, float znear, float zfar);
    // Constructs a perspective projection matrix (compatible with OpenGL)
    static Mat4 Frustum(float left, float right, float bottom, float top, float znear, float zfar);
    // Constructs a perspective projection matrix (compatible with OpenGL)
    static Mat4 Perspective(float fov, float aspect, float znear, float zfar);
};

// Quaternion representation of rotations.
// R(v) = M v = q v q'
struct Quat
{
    float w; // scalar part
    float x; // (bi-)vector part
    float y;
    float z;

    Quat() = default;

    Quat(float w, float x, float y, float z);

    Quat(float w, Vec3 const& v);

    // Construct from a 3D vector, i.e. construct a pure quaternion.
    explicit Quat(Vec3 const& v);

    // Returns the vector part
    explicit operator Vec3() const;

    // Returns the scalar part
    float scalar() const;

    // Returns the vector part
    Vec3 vec() const;

    // Returns the unit quaternion
    static Quat Identity();
    // Constructs a quaternion, which represents a rotation around the given axis
    static Quat Rotation(Vec3 const& axis, float angle);
    // Constructs a quaternion, which rotates FROM into TO.
    static Quat Rotation(Vec3 const& from, Vec3 const& to);
    // Construct a quaternion from Euler angles.
    static Quat Rotation(Euler const& e);
    // Constructs a quaternion q from a matrix M such that q v q' = M v.
    static Quat Rotation(Mat3 const& m);
    // Constructs a quaternion q from a matrix M such that q v q' = M v.
    // NOTE: Only the rotational part of the matrix is used.
    static Quat Rotation(Mat3x4 const& m);
    // Constructs a quaternion q from a matrix M such that q v q' = M v.
    // NOTE: Only the rotational part of the matrix is used.
    static Quat Rotation(Mat4 const& m);
};

// Euler angle representation of rotations.
// Order is:
//  1. heading (y-axis)
//  2. elevation (x-axis)
//  3. bank (z-axis)
struct Euler
{
    float heading;
    float elevation;
    float bank;

    Euler() = default;

    Euler(float heading, float elevation, float bank);

    // Identity transform.
    static Euler Identity();
    // Extract Euler angles from a rotation matrix
    static Euler Rotation(Mat3 const& m);
    // Extract Euler angles from a quaternion
    static Euler Rotation(Quat const& q);
};

//==============================================================================
// Matrix arithmetic
//
// Note:
// If the matrices don't match, use Mul()
//==============================================================================

// Returns the cross-product matrix for u, i.e. returns a matrix M_u, such that
// M_u * v = Cross(u,v) for all v.
inline Mat3 CrossProductMatrix(Vec3 const& u) {
    return { 0, u.z, -u.y, -u.z, 0, u.x, u.y, -u.x, 0 };
}

// Returns u v^T
inline Mat3 OuterProduct(Vec3 const& u, Vec3 const& v) {
    return { u * v.x, u * v.y, u * v.z };
}

// Returns u v^T
inline Mat4 OuterProduct(Vec4 const& u, Vec4 const& v) {
    return { u * v.x, u * v.y, u * v.z, u * v.w };
}

//------------------------------------------------------------------------------
// Matrix * Scalar
//

inline Mat3   operator *(Mat3   const& A, float s) { return { A.col0 * s, A.col1 * s, A.col2 * s }; }
inline Mat3x4 operator *(Mat3x4 const& A, float s) { return { A.col0 * s, A.col1 * s, A.col2 * s, A.col3 * s }; }
inline Mat4   operator *(Mat4   const& A, float s) { return { A.col0 * s, A.col1 * s, A.col2 * s, A.col3 * s }; }

inline Mat3   operator *(float s, Mat3   const& A) { return A * s; }
inline Mat3x4 operator *(float s, Mat3x4 const& A) { return A * s; }
inline Mat4   operator *(float s, Mat4   const& A) { return A * s; }

//------------------------------------------------------------------------------
// Matrix * Vector
//

template <class Mat, class Vec> auto _lin2(Mat const& M, Vec const& v) { return M.col0 * v.x + M.col1 * v.y; }
template <class Mat, class Vec> auto _lin3(Mat const& M, Vec const& v) { return M.col0 * v.x + M.col1 * v.y + M.col2 * v.z; }
template <class Mat, class Vec> auto _lin4(Mat const& M, Vec const& v) { return M.col0 * v.x + M.col1 * v.y + M.col2 * v.z + M.col3 * v.w; }

// Returns A * v
inline Vec3 operator *(Mat3 const& A, Vec3 const& v) { return _lin3(A, v); }

// Returns A * v
inline Point3 operator *(Mat3 const& A, Point3 const& v)
{
    // [A 0] [v] = [A*v]
    // [0 1] [1]   [1  ]

    return Point3(_lin3(A, v));
}

// Returns Mat4(A) * v
inline Vec4 Mul(Mat3 const& A, Vec4 const& v)
{
    // [A 0] [v] = [A*v]
    // [0 1] [s]   [s  ]

    return Vec4(_lin3(A, v), v.w);
}

// Returns A * (x,y,z,0)
inline Vec3 Mul(Mat3x4 const& A, Vec3 const& v)
{
    // [A a] [v] = [A*v]
    // [0 1] [0]   [1  ]

    return _lin3(A, v);
}

// Returns A * (x,y,z,1)
inline Point3 Mul(Mat3x4 const& A, Point3 const& p)
{
    // [A a] [v] = [A*v+a]
    // [0 1] [1]   [1    ]

    return Point3(_lin3(A, p) + A.col3);
}

// Returns A * v
inline Vec3 operator *(Mat3x4 const& A, Vec4 const& v) { return _lin4(A, v); }

// Returns A * (x,y,z,0)
inline Vec4 Mul(Mat4 const& A, Vec3 const& v)
{
    // [A a] [v] = [A*v]
    // [b c] [0]   [b*v]

    return _lin3(A, v);
}

// Returns A * (x,y,z,1) projected back onto w = 1.
inline Point3 Mul(Mat4 const& A, Point3 const& p)
{
    // [A a] [v] = [A*v+a]
    // [b c] [1]   [b*v+c]

    return Point3(_lin3(A, p) + A.col3);
}

// Returns A * v
inline Vec4 operator *(Mat4 const& A, Vec4 const& v) { return _lin4(A, v); }

//------------------------------------------------------------------------------
// Matrix * Matrix
//

// Returns A * B
inline Mat3 operator *(Mat3 const& A, Mat3 const& B) {
    return { A * B.col0, A * B.col1, A * B.col2 };
}

// Returns A * B
inline Mat3x4 operator *(Mat3 const& A, Mat3x4 const& B) {
    return { A * B.col0, A * B.col1, A * B.col2, A * B.col3 };
}

// Returns Mat4(A) * B
inline Mat4 Mul(Mat3 const& A, Mat4 const& B)
{
    // [A 0] [B b] = [A*B A*b]
    // [0 1] [c d]   [c   d  ]

    const auto c0 = Vec4(_lin3(A, B.col0), B.col0.w);
    const auto c1 = Vec4(_lin3(A, B.col1), B.col1.w);
    const auto c2 = Vec4(_lin3(A, B.col2), B.col2.w);
    const auto c3 = Vec4(_lin3(A, B.col3), B.col3.w);

    return { c0, c1, c2, c3 };
}

// Returns A * Mat4(B)
inline Mat3x4 Mul(Mat3x4 const& A, Mat3 const& B)
{
    // [A a] [B 0] = [A*B a]
    // [0 1] [0 1]   [0   1]

    const auto c0 = _lin3(A, B.col0);
    const auto c1 = _lin3(A, B.col1);
    const auto c2 = _lin3(A, B.col2);

    return { c0, c1, c2, A.col3 };
}

// Returns A * Mat4(B)
inline Mat3x4 Mul(Mat3x4 const& A, Mat3x4 const& B)
{
    // [A a] [B b] = [A*B A*b+a]
    // [0 1] [0 1]   [0   1    ]

    const auto c0 = _lin3(A, B.col0);
    const auto c1 = _lin3(A, B.col1);
    const auto c2 = _lin3(A, B.col2);
    const auto c3 = _lin3(A, B.col3) + A.col3;

    return { c0, c1, c2, c3 };
}

// Returns A * B
inline Mat3x4 operator *(Mat3x4 const& A, Mat4 const& B) {
    return { A * B.col0, A * B.col1, A * B.col2, A * B.col3 };
}

// Return A * Mat4(B)
inline Mat4 Mul(Mat4 const& A, Mat3 const& B)
{
    // [A a] [B 0] = [A*B a]
    // [u v] [0 1]   [u*B v]

    const auto c0 = _lin3(A, B.col0);
    const auto c1 = _lin3(A, B.col1);
    const auto c2 = _lin3(A, B.col2);

    return { c0, c1, c2, A.col3 };
}

// Return A * Mat4(B)
inline Mat4 Mul(Mat4 const& A, Mat3x4 const& B)
{
    // [A a] [B b] = [A*B A*b+a]
    // [u v] [0 1]   [u*B u*b+v]

    const auto c0 = _lin3(A, B.col0);
    const auto c1 = _lin3(A, B.col1);
    const auto c2 = _lin3(A, B.col2);
    const auto c3 = _lin3(A, B.col3) + A.col3;

    return { c0, c1, c2, c3 };
}

// Return A * B
inline Mat4 operator *(Mat4 const& A, Mat4 const& B) {
    return { A * B.col0, A * B.col1, A * B.col2, A * B.col3 };
}

//==============================================================================
// Mat3
//==============================================================================

inline Mat3::Mat3(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2)
    : col0(c0)
    , col1(c1)
    , col2(c2)
{
}

inline Mat3::Mat3(
        float m00, float m10, float m20,
        float m01, float m11, float m21,
        float m02, float m12, float m22)
    : col0(m00, m10, m20)
    , col1(m01, m11, m21)
    , col2(m02, m12, m22)
{
}

inline Mat3::Mat3(float const* p)
    : col0(p + 0)
    , col1(p + 3)
    , col2(p + 6)
{
}

inline Mat3::Mat3(float m00, float m11, float m22)
    : col0(m00, 0.0f, 0.0f)
    , col1(0.0f, m11, 0.0f)
    , col2(0.0f, 0.0f, m22)
{
}

inline Mat3::Mat3(Mat3x4 const& A)
    : col0(A.col0)
    , col1(A.col1)
    , col2(A.col2)
{
}

inline Mat3::Mat3(Mat4 const& A)
    : col0(A.col0)
    , col1(A.col1)
    , col2(A.col2)
{
}

inline float* Mat3::data() {
    return reinterpret_cast<float*>(this);
}

inline float const* Mat3::data() const {
    return reinterpret_cast<float const*>(this);
}

inline Vec3& Mat3::operator()(int col) {
    return (&col0)[col];
}

inline Vec3 const& Mat3::operator()(int col) const {
    return (&col0)[col];
}

inline float& Mat3::operator()(int row, int col) {
    return (operator()(col))[row];
}

inline float const& Mat3::operator()(int row, int col) const {
    return (operator()(col))[row];
}

inline Mat3 Mat3::FromRows(Vec3 const& r0, Vec3 const& r1, Vec3 const& r2) {
    return { r0.x, r1.x, r2.x, r0.y, r1.y, r2.y, r0.z, r1.z, r2.z };
}

inline Mat3 Mat3::Identity() {
    return { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
}

inline Mat3 Mat3::Zero() {
    return { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
}

inline Mat3 Mat3::Shear(float xy, float yx, float xz, float zx, float yz, float zy) {
    return { 1, xy, xz, yx, 1, yz, zx, zy, 1};
}

inline Mat3 Mat3::Rotation(Vec3 const& axis, float angle)
{
    const auto n = Normalize(axis);
    const auto s = Sin(angle);
    const auto c = Cos(angle);

    const auto k = 1.0f - c;

    return {
        n.x * n.x * k + c,
        n.x * n.y * k + s * n.z,
        n.x * n.z * k - s * n.y,
        n.y * n.x * k - s * n.z,
        n.y * n.y * k + c,
        n.y * n.z * k + s * n.x,
        n.z * n.x * k + s * n.y,
        n.z * n.y * k - s * n.x,
        n.z * n.z * k + c
        };
}

inline Mat3 Mat3::Rotation(Vec3 const& from, Vec3 const& to)
{
    //
    // XXX:
    // Use Quat::rotation(from, to) ???
    //

    const auto f = Normalize(from);
    const auto t = Normalize(to);

    const auto n = Cross(f, t);
    const auto e = Dot(f, t);
    const auto h = 1.0f / (1.0f + e);

    return {
        h * n.x * n.x + e,
        h * n.x * n.y + n.z,
        h * n.x * n.z - n.y,
        h * n.x * n.y - n.z,
        h * n.y * n.y + e,
        h * n.y * n.z + n.x,
        h * n.x * n.z + n.y,
        h * n.y * n.z - n.x,
        h * n.z * n.z + e
        };
}

inline Mat3 Mat3::Rotation(Quat const& q)
{
    return {
        2.0f * (q.w * q.w + q.x * q.x) - 1.0f,
        2.0f * (q.x * q.y + q.w * q.z),
        2.0f * (q.x * q.z - q.w * q.y),
        2.0f * (q.x * q.y - q.w * q.z),
        2.0f * (q.w * q.w + q.y * q.y) - 1.0f,
        2.0f * (q.y * q.z + q.w * q.x),
        2.0f * (q.x * q.z + q.w * q.y),
        2.0f * (q.y * q.z - q.w * q.x),
        2.0f * (q.w * q.w + q.z * q.z) - 1.0f
        };
}

inline Mat3 Mat3::Rotation(Euler const& e) {
    return Mat3::Rotation(Quat::Rotation(e));
}

inline Mat3 Mat3::Frame(Vec3 const& n)
{
    const auto dz = Normalize(n);
    const auto dx = Normalize(Perp(dz));
    const auto dy = Cross(dz, dx); // normalize?!?!

    return { dx, dy, dz };
}

inline Vec3 Mat3::ToLocal(Vec3 const& v) const {
    return { Dot(col0, v), Dot(col1, v), Dot(col2, v) };
}

inline Vec3 Mat3::ToWorld(Vec3 const& v) const {
    return col0 * v.x + col1 * v.y + col2 * v.z;
}

inline float Mat3::CosTheta2(Vec3 const& v)
{
//  return Saturate(v.z * v.z);
    return v.z * v.z;
}

inline float Mat3::CosTheta(Vec3 const& v)
{
//  return Clamp(v.z, -1.0f, 1.0f);
    return v.z;
}

inline float Mat3::SinTheta2(Vec3 const& v)
{
//  return Saturate(1.0f - v.z * v.z);
    return 1.0f - v.z * v.z;
}

inline float Mat3::SinTheta(Vec3 const& v)
{
//  return Sqrt(SinTheta2(v));
    return Sqrt(Max0(1.0f - v.z * v.z));
}

inline float Mat3::TanTheta2(Vec3 const& v)
{
//  return SinTheta2(v) / (v.z * v.z);
    return (1.0f - v.z * v.z) / (v.z * v.z);
}

inline float Mat3::TanTheta(Vec3 const& v)
{
//  return SinTheta(v) / v.z;
    return Sqrt(Max0(1.0f - v.z * v.z)) / v.z;
}

inline float Mat3::SinPhi2(Vec3 const& v)
{
//  return Saturate(v.y * v.y / SinTheta2(v));
    return v.y * v.y / SinTheta2(v);
}

inline float Mat3::SinPhi(Vec3 const& v)
{
//  return Clamp(v.y / SinTheta(v), -1.0f, 1.0f);
    return v.y / SinTheta(v);
}

inline float Mat3::CosPhi2(Vec3 const& v)
{
//  return Saturate(v.x * v.x / SinTheta2(v));
    return v.x * v.x / SinTheta2(v);
}

inline float Mat3::CosPhi(Vec3 const& v)
{
//  return Clamp(v.x / SinTheta(v), -1.0f, 1.0f);
    return v.x / SinTheta(v);
}

//------------------------------------------------------------------------------
// Geometric functions
//

inline Mat3 Transpose(Mat3 const& m) {
    return Mat3::FromRows(m.col0, m.col1, m.col2);
}

template <class Mat>
inline Mat3 _inv3(Mat const& A)
{
    const auto r0 = Cross(A.col1, A.col2);
    const auto r1 = Cross(A.col2, A.col0);
    const auto r2 = Cross(A.col0, A.col1);

    const auto det = Dot(A.col0, r0);

    return Mat3::FromRows(r0 / det, r1 / det, r2 / det);
}

// Returns a linear transform B such that B(A(x)) = x
inline Mat3 Inverse(Mat3 const& A) {
    return _inv3(A);
}

// Returns the determinant of the matrix
inline float Determinant(Mat3 const& A) {
    return Dot(A.col0, Cross(A.col1, A.col2));
}

// Constructs a linear map R from a normal N such that R*(0,0,1) = N.
//inline Mat3 Frame(Vec3 const& n) {
//    return Mat3::frame(n);
//}

template <class Mat>
float _rangle(Mat const& A) // Rotation angle.
{
    return SafeAcos(0.5f * (A(0,0) + A(1,1) + A(2,2) - 1.0f));
}

template <class Mat>
Vec3 _raxis(Mat const& A) // Rotation axis. NOT normalized.
{
    return Vec3(A(2,1) - A(1,2), A(0,2) - A(2,0), A(1,0) - A(0,1));
}

// Returns the rotation angle for the given rotation matrix.
inline float RotationAngle(Mat3 const& A) {
    return _rangle(A);
}

// Returns the rotation axis for the rotation matrix.
inline Vec3 RotationAxis(Mat3 const& A) {
    return _raxis(A);
}

// Factor A = Q * R using modified Gram-Schmidt, where
// Q is an orthonormal matrix and R and upper-triangular matrix.
//
//struct QR {
//    Mat3 Q;
//    Mat3 R; // upper triangular
//};
//
//inline QR FactorQR(Mat3 const& A)
//{
//    QR qr;
//
//    auto v0 = A.col0;
//    auto v1 = A.col1;
//    auto v2 = A.col2;
//
//    auto r00  = Length(v0);
//    qr.Q.col0 = v0 / r00;
//    auto r01  = Dot(qr.Q.col0, v1);
//    auto r02  = Dot(qr.Q.col0, v2);
//    v1        = v1 - r01 * qr.Q.col0;
//    v2        = v2 - r02 * qr.Q.col0;
//
//    auto r11  = Length(v1);
//    qr.Q.col1 = v1 / r11;
//    auto r12  = Dot(qr.Q.col1, v2);
//    v2        = v2 - r12 * qr.Q.col1;
//
//    auto r22  = Length(v2);
//    qr.Q.col2 = v2 / r22;
//
//    qr.R.col0 = {r00,   0,   0};
//    qr.R.col1 = {r01, r11,   0};
//    qr.R.col2 = {r02, r12, r22};
//
//    return qr;
//}

// Orthonormalize the given matrix A.
//inline Mat3 Orthonormalize(Mat3 const& A)
//{
//#if 0
//    QR qr = FactorQR(A);
//    return qr.Q;
//#else
//    Mat3 Q;
//
//    Q.col0 = Normalize(A.col0);
//    Q.col1 = Normalize(A.col1 - Dot(Q.col0, A.col1) * Q.col0);
//
//#if 0
//    // This is valid for right-handed coordinate systems only...
//    Q.col2 = Cross(Q.col0, Q.col1);
//#else
//    auto v = A.col2 - Dot(Q.col0, A.col2) * Q.col0;
//    Q.col2 = Normalize(v - Dot(Q.col1, v) * Q.col1);
//#endif
//
//    return Q;
//#endif
//}

//==============================================================================
// Mat3x4
//==============================================================================

inline Mat3x4::Mat3x4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2, Vec3 const& c3)
    : col0(c0)
    , col1(c1)
    , col2(c2)
    , col3(c3)
{
}

inline Mat3x4::Mat3x4(
        float m00, float m10, float m20,
        float m01, float m11, float m21,
        float m02, float m12, float m22,
        float m03, float m13, float m23)
    : col0(m00, m10, m20)
    , col1(m01, m11, m21)
    , col2(m02, m12, m22)
    , col3(m03, m13, m23)
{
}

inline Mat3x4::Mat3x4(Mat3 const& L)
    : col0(L.col0)
    , col1(L.col1)
    , col2(L.col2)
    , col3(0.0f)
{
}

inline Mat3x4::Mat3x4(Mat3 const& L, Vec3 const& T)
    : col0(L.col0)
    , col1(L.col1)
    , col2(L.col2)
    , col3(T)
{
}

inline Mat3x4::Mat3x4(Mat3 const& L, Point3 const& T)
    : col0(L.col0)
    , col1(L.col1)
    , col2(L.col2)
    , col3(T.x, T.y, T.z)
{
}

inline Mat3x4::Mat3x4(Mat4 const& A)
    : col0(A.col0)
    , col1(A.col1)
    , col2(A.col2)
    , col3(A.col3)
{
}

inline float* Mat3x4::data() {
    return reinterpret_cast<float*>(this);
}

inline float const* Mat3x4::data() const {
    return reinterpret_cast<float const*>(this);
}

inline Vec3& Mat3x4::operator()(int col) {
    return (&col0)[col];
}

inline Vec3 const& Mat3x4::operator()(int col) const {
    return (&col0)[col];
}

inline float& Mat3x4::operator()(int row, int col) {
    return (operator()(col))[row];
}

inline float const& Mat3x4::operator()(int row, int col) const {
    return (operator()(col))[row];
}

inline Mat3x4 Mat3x4::Identity() {
    return Mat3x4(Mat3::Identity());
}

inline Mat3x4 Mat3x4::Zero() {
    return Mat3x4(Mat3::Zero());
}

inline Mat3x4 Mat3x4::Translation(Vec3 const& t) {
    return Mat3x4(Mat3::Identity(), t);
}

inline Mat3x4 Mat3x4::Scaling(Vec3 const& s) {
    return Mat3x4(Mat3::Scaling(s));
}

inline Mat3x4 Mat3x4::Scaling(Vec3 const& s, Vec3 const& origin)
{
    Mat3x4 out;

    out.col0 = {s.x, 0, 0};
    out.col1 = {0, s.y, 0};
    out.col2 = {0, 0, s.z};
    out.col3 = origin * (1.0f - s);

    return out;
}

inline Mat3x4 Mat3x4::Shear(float xy, float yx, float xz, float zx, float yz, float zy) {
    return Mat3x4(Mat3::Shear(xy, yx, xz, zx, yz, zy));
}

inline Mat3x4 Mat3x4::Rotation(Vec3 const& axis, float angle) {
    return Mat3x4(Mat3::Rotation(axis, angle));
}

inline Mat3x4 Mat3x4::Rotation(Vec3 const& from, Vec3 const& to) {
    return Mat3x4(Mat3::Rotation(from, to));
}

inline Mat3x4 Mat3x4::Rotation(Quat const& q) {
    return Mat3x4(Mat3::Rotation(q));
}

inline Mat3x4 Mat3x4::Rotation(Euler const& e) {
    return Mat3x4(Mat3::Rotation(e));
}

inline Mat3x4 Mat3x4::LookAt(Point3 const& position, Point3 const& target, Vec3 const& up)
{
    const auto z = Normalize(position - target);
    const auto x = Normalize(Cross(up, z));
    const auto y = Cross(z, x); // normalize?!?!

    return {
        x.x,
        y.x,
        z.x,
        x.y,
        y.y,
        z.y,
        x.z,
        y.z,
        z.z,
       -Dot(position, x),
       -Dot(position, y),
       -Dot(position, z)
        };
}

//------------------------------------------------------------------------------
// Geometric functions
//

// Returns an affine transform T such that T(A(x)) = x
inline Mat3x4 Inverse(Mat3x4 const& A)
{
    const auto Linv = _inv3(A);
    return Mat3x4(Linv, -(Linv * A.col3));
}

// Returns the rotation angle for the given rotation matrix.
// Only the upper-left 3x3 submatrix is used and assumed to be a proper rotation matrix.
inline float RotationAngle(Mat3x4 const& A) {
    return _rangle(A);
}

// Returns the rotation axis for the rotation matrix.
// Only the upper-left 3x3 submatrix is used and assumed to be a proper rotation matrix.
inline Vec3 RotationAxis(Mat3x4 const& A) {
    return _raxis(A);
}

//==============================================================================
// Mat4
//==============================================================================

inline Mat4::Mat4(Vec4 const& c0, Vec4 const& c1, Vec4 const& c2, Vec4 const& c3)
    : col0(c0)
    , col1(c1)
    , col2(c2)
    , col3(c3)
{
}

inline Mat4::Mat4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2)
    : col0(c0, 0.0f)
    , col1(c1, 0.0f)
    , col2(c2, 0.0f)
    , col3(0.0f, 0.0f, 0.0f, 1.0f)
{
}

inline Mat4::Mat4(Vec3 const& c0, Vec3 const& c1, Vec3 const& c2, Vec3 const& c3)
    : col0(c0, 0.0f)
    , col1(c1, 0.0f)
    , col2(c2, 0.0f)
    , col3(c3, 1.0f)
{
}

inline Mat4::Mat4(
        float m00, float m10, float m20, float m30,
        float m01, float m11, float m21, float m31,
        float m02, float m12, float m22, float m32,
        float m03, float m13, float m23, float m33)
    : col0(m00, m10, m20, m30)
    , col1(m01, m11, m21, m31)
    , col2(m02, m12, m22, m32)
    , col3(m03, m13, m23, m33)
{
}

inline Mat4::Mat4(float m00, float m11, float m22, float m33)
    : col0(m00, 0.0f, 0.0f, 0.0f)
    , col1(0.0f, m11, 0.0f, 0.0f)
    , col2(0.0f, 0.0f, m22, 0.0f)
    , col3(0.0f, 0.0f, 0.0f, m33)
{
}

inline Mat4::Mat4(float const* p)
    : col0(p +  0)
    , col1(p +  4)
    , col2(p +  8)
    , col3(p + 12)
{
}

inline Mat4::Mat4(Mat3 const& L)
    : col0(L.col0, 0.0f)
    , col1(L.col1, 0.0f)
    , col2(L.col2, 0.0f)
    , col3(0.0f, 0.0f, 0.0f, 1.0f)
{
}

inline Mat4::Mat4(Mat3 const& L, Vec3 const& T)
    : col0(L.col0, 0.0f)
    , col1(L.col1, 0.0f)
    , col2(L.col2, 0.0f)
    , col3(T.x, T.y, T.z, 1.0f)
{
}

inline Mat4::Mat4(Mat3 const& L, Point3 const& T)
    : col0(L.col0, 0.0f)
    , col1(L.col1, 0.0f)
    , col2(L.col2, 0.0f)
    , col3(T.x, T.y, T.z, 1.0f)
{
}

inline Mat4::Mat4(Mat3x4 const& A)
    : col0(A.col0, 0.0f)
    , col1(A.col1, 0.0f)
    , col2(A.col2, 0.0f)
    , col3(A.col3, 1.0f)
{
}

inline float* Mat4::data() {
    return reinterpret_cast<float*>(this);
}

inline float const* Mat4::data() const {
    return reinterpret_cast<float const*>(this);
}

inline Vec4& Mat4::operator()(int col) {
    return (&col0)[col];
}

inline Vec4 const& Mat4::operator()(int col) const {
    return (&col0)[col];
}

inline float& Mat4::operator()(int row, int col) {
    return (operator()(col))[row];
}

inline float const& Mat4::operator()(int row, int col) const {
    return (operator()(col))[row];
}

inline Mat4 Mat4::FromRows(Vec4 const& r0, Vec4 const& r1, Vec4 const& r2, Vec4 const& r3)
{
    return {
        r0.x, r1.x, r2.x, r3.x,
        r0.y, r1.y, r2.y, r3.y,
        r0.z, r1.z, r2.z, r3.z,
        r0.w, r1.w, r2.w, r3.w
        };
}

inline Mat4 Mat4::Identity() {
    return { 1.0f, 1.0f, 1.0f, 1.0f };
}

inline Mat4 Mat4::Zero() {
    return { 0.0f, 0.0f, 0.0f, 0.0f };
}

inline Mat4 Mat4::Translation(Vec3 const& t) {
    return Mat4(Mat3x4::Translation(t));
}

inline Mat4 Mat4::Scaling(Vec3 const& s) {
    return Mat4(Mat3::Scaling(s));
}

inline Mat4 Mat4::Scaling(Vec3 const& s, Vec3 const& origin)
{
    Mat4 out;

    out.col0 = {s.x, 0, 0, 0};
    out.col1 = {0, s.y, 0, 0};
    out.col2 = {0, 0, s.z, 0};
    out.col3 = Vec4(origin * (1.0f - s), 1);

    return out;
}

inline Mat4 Mat4::Shear(float xy, float yx, float xz, float zx, float yz, float zy) {
    return Mat4(Mat3::Shear(xy, yx, xz, zx, yz, zy));
}

inline Mat4 Mat4::Rotation(Vec3 const& axis, float angle) {
    return Mat4(Mat3::Rotation(axis, angle));
}

inline Mat4 Mat4::Rotation(Vec3 const& from, Vec3 const& to) {
    return Mat4(Mat3::Rotation(from, to));
}

inline Mat4 Mat4::Rotation(Quat const& q) {
    return Mat4(Mat3::Rotation(q));
}

inline Mat4 Mat4::Rotation(Euler const& e) {
    return Mat4(Mat3::Rotation(e));
}

inline Mat4 Mat4::LookAt(Point3 const& position, Point3 const& target, Vec3 const& up) {
    return Mat4(Mat3x4::LookAt(position, target, up));
}

inline Mat4 Mat4::Ortho(float left, float right, float bottom, float top, float znear, float zfar)
{
    const auto idx = 1.0f / (right - left);
    const auto idy = 1.0f / (top - bottom);
    const auto idz = 1.0f / (zfar - znear);

    return {
        2 * idx,
        0,
        0,
        0,
        0,
        2 * idy,
        0,
        0,
        0,
        0,
        2 * idz,
        0,
       -(right + left) * idx,
       -(top + bottom) * idy,
       -(zfar + znear) * idz,
        1
        };
}

inline Mat4 Mat4::Frustum(float left, float right, float bottom, float top, float znear, float zfar)
{
    const auto idx = 1.0f / (right - left);
    const auto idy = 1.0f / (top - bottom);
    const auto idz = 1.0f / (zfar - znear);

    return {
        2 * znear * idx,
        0,
        0,
        0,
        0,
        2 * znear * idy,
        0,
        0,
        (right + left) * idx,
        (top + bottom) * idy,
       -(zfar + znear) * idz,
       -1,
        0,
        0,
       -2 * zfar * znear * idz,
        0
        };
}

inline Mat4 Mat4::Perspective(float fov, float aspect, float znear, float zfar)
{
    const auto f = 1.0f / Tan(0.5f * fov);

    return {
        f / aspect,
        0,
        0,
        0,
        0,
        f,
        0,
        0,
        0,
        0,
       -(zfar + znear) / (zfar - znear),
       -1,
        0,
        0,
       -2 * zfar * znear / (zfar - znear),
        0
        };
}

//------------------------------------------------------------------------------
// Geometric functions
//

// Returns the transpose of the given matrix A
inline Mat4 Transpose(Mat4 const& A) {
    return Mat4::FromRows(A.col0, A.col1, A.col2, A.col3);
}

// Returns the inverse of the given matrix A.
inline Mat4 Inverse(Mat4 const& A)
{
    const auto s0 = Det2(A(0,0), A(0,1), A(1,0), A(1,1));
    const auto s1 = Det2(A(0,0), A(0,2), A(1,0), A(1,2));
    const auto s2 = Det2(A(0,0), A(0,3), A(1,0), A(1,3));
    const auto s3 = Det2(A(0,1), A(0,2), A(1,1), A(1,2));
    const auto s4 = Det2(A(0,1), A(0,3), A(1,1), A(1,3));
    const auto s5 = Det2(A(0,2), A(0,3), A(1,2), A(1,3));
    const auto c5 = Det2(A(2,2), A(2,3), A(3,2), A(3,3));
    const auto c4 = Det2(A(2,1), A(2,3), A(3,1), A(3,3));
    const auto c3 = Det2(A(2,1), A(2,2), A(3,1), A(3,2));
    const auto c2 = Det2(A(2,0), A(2,3), A(3,0), A(3,3));
    const auto c1 = Det2(A(2,0), A(2,2), A(3,0), A(3,2));
    const auto c0 = Det2(A(2,0), A(2,1), A(3,0), A(3,1));

    const auto idet = 1.0f / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

    return {
        idet * (+A(1,1) * c5 - A(1,2) * c4 + A(1,3) * c3),
        idet * (-A(1,0) * c5 + A(1,2) * c2 + A(1,3) * c1),
        idet * (+A(1,0) * c4 - A(1,1) * c2 + A(1,3) * c0),
        idet * (-A(1,0) * c3 + A(1,1) * c1 + A(1,2) * c0),
        idet * (-A(0,1) * c5 + A(0,2) * c4 - A(0,3) * c3),
        idet * (+A(0,0) * c5 - A(0,2) * c2 + A(0,3) * c1),
        idet * (-A(0,0) * c4 + A(0,1) * c2 - A(0,3) * c0),
        idet * (+A(0,0) * c3 - A(0,1) * c1 + A(0,2) * c0),
        idet * (+A(3,1) * s5 - A(3,2) * s4 + A(3,3) * s3),
        idet * (-A(3,0) * s5 + A(3,2) * s2 - A(3,3) * s1),
        idet * (+A(3,0) * s4 - A(3,1) * s2 + A(3,3) * s0),
        idet * (-A(3,0) * s3 + A(3,1) * s1 - A(3,2) * s0),
        idet * (-A(2,1) * s5 + A(2,2) * s4 - A(2,3) * s3),
        idet * (+A(2,0) * s5 - A(2,2) * s2 + A(2,3) * s1),
        idet * (-A(2,0) * s4 + A(2,1) * s2 - A(2,3) * s0),
        idet * (+A(2,0) * s3 - A(2,1) * s1 + A(2,2) * s0)
        };
}

// Returns the rotation angle for the given rotation matrix.
// Only the upper-left 3x3 submatrix is used and assumed to be a proper rotation matrix.
inline float RotationAngle(Mat4 const& A) {
    return _rangle(A);
}

// Returns the rotation axis for the rotation matrix.
// Only the upper-left 3x3 submatrix is used and assumed to be a proper rotation matrix.
inline Vec3 RotationAxis(Mat4 const& A) {
    return _raxis(A);
}

//==============================================================================
// Quat
//==============================================================================

inline Quat::Quat(float w, float x, float y, float z)
    : w(w)
    , x(x)
    , y(y)
    , z(z)
{
}

inline Quat::Quat(float w, Vec3 const& v)
    : w(w)
    , x(v.x)
    , y(v.y)
    , z(v.z)
{
}

inline Quat::Quat(Vec3 const& v)
    : w(0)
    , x(v.x)
    , y(v.y)
    , z(v.z)
{
}

inline Quat::operator Vec3() const {
    return {x, y, z};
}

inline float Quat::scalar() const {
    return w;
}

inline Vec3 Quat::vec() const {
    return {x, y, z};
}

inline Quat Quat::Identity() {
    return { 1, 0, 0, 0 };
}

inline Quat Quat::Rotation(Vec3 const& axis, float angle)
{
    float s, c;
    SinCos(0.5f * angle, s, c);

    s /= Length(axis);

    return { c, s * axis.x, s * axis.y, s * axis.z };
}

inline Quat Quat::Rotation(Vec3 const& from, Vec3 const& to)
{
    // Result is normalized if FROM and TO are normalized.

    const auto d = Clamp(Dot(from, to), -1.0f, 1.0f);

#if 0
    if (d == -1.0f)
        return Quat(0.0f, Normalize(from)); // Quat::Rotation(Perp(from), kPi)
#endif

    const auto f = Sqrt(2.0f + 2.0f * d);
    const auto v = Cross(from, to);

    return { 0.5f * f, v.x / f, v.y / f, v.z / f };
}

inline Quat Quat::Rotation(Euler const& e)
{
    // 1: heading (y)
    // 2: elevation (x)
    // 3: bank (z)

    float sh, ch, se, ce, sb, cb;
    SinCos(0.5f * e.heading,   sh, ch);
    SinCos(0.5f * e.elevation, se, ce);
    SinCos(0.5f * e.bank,      sb, cb);

    const float w = ch * ce * cb - sh * se * sb;
    const float x = ch * se * cb - sh * ce * sb;
    const float y = sh * ce * cb + ch * se * sb;
    const float z = ch * ce * sb + sh * se * cb;

    return Quat(w, x, y, z);
}

template <class Mat>
Quat _mtoq(Mat const& A)
{
    float w = 0.5f * SafeSqrt(1.0f + A(0,0) + A(1,1) + A(2,2));
    float x = 0.5f * SafeSqrt(1.0f + A(0,0) - A(1,1) - A(2,2));
    float y = 0.5f * SafeSqrt(1.0f - A(0,0) + A(1,1) - A(2,2));
    float z = 0.5f * SafeSqrt(1.0f - A(0,0) - A(1,1) + A(2,2));

    float m = Max(w, x, y, z);
    if (m == w)
    {
        x = Copysign(x, A(2,1) - A(1,2));
        y = Copysign(y, A(0,2) - A(2,0));
        z = Copysign(z, A(1,0) - A(0,1));
    }
    else if (m == x)
    {
        w = Copysign(w, A(2,1) - A(1,2));
        y = Copysign(y, A(1,0) + A(0,1));
        z = Copysign(z, A(0,2) + A(2,0));
    }
    else if (m == y)
    {
        w = Copysign(w, A(0,2) - A(2,0));
        x = Copysign(x, A(1,0) + A(0,1));
        z = Copysign(z, A(2,1) + A(1,2));
    }
    else // m == z
    {
        w = Copysign(w, A(1,0) - A(0,1));
        x = Copysign(x, A(2,0) + A(0,2));
        y = Copysign(y, A(2,1) + A(1,2));
    }

#if 0
    return Normalize(Quat(w, x, y, z));
#else
    const float len = Hypot(w, x, y, z);
    return { w / len, x / len, y / len, z / len };
#endif
}

inline Quat Quat::Rotation(Mat3   const& A) { return _mtoq(A); }
inline Quat Quat::Rotation(Mat3x4 const& A) { return _mtoq(A); }
inline Quat Quat::Rotation(Mat4   const& A) { return _mtoq(A); }

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator==(Quat const& p, Quat const& q) {
    return p.w == q.w && p.x == q.x && p.y == q.y && p.z == q.z;
}

inline bool operator!=(Quat const& p, Quat const& q) {
    return !(p == q);
}

//------------------------------------------------------------------------------
// Arithmetic
//

inline Quat operator -(Quat const& p) {
    return { -p.w, -p.x, -p.y, -p.z };
}

inline Quat operator +(Quat const& p, Quat const& q) {
    return { p.w + q.w, p.x + q.x, p.y + q.y, p.z + q.z };
}

inline Quat operator -(Quat const& p, Quat const& q) {
    return { p.w - q.w, p.x - q.x, p.y - q.y, p.z - q.z };
}

inline Quat operator *(Quat const& p, Quat const& q)
{
    const float w = p.w * q.w - p.x * q.x - p.y * q.y - p.z * q.z;
    const float x = p.w * q.x + p.x * q.w + p.y * q.z - p.z * q.y;
    const float y = p.w * q.y + p.y * q.w + p.z * q.x - p.x * q.z;
    const float z = p.w * q.z + p.z * q.w + p.x * q.y - p.y * q.x;

    return Quat(w, x, y, z);
}

inline Quat operator *(Quat const& p, float s) {
    return { p.w * s, p.x * s, p.y * s, p.z * s };
}

inline Quat operator /(Quat const& p, float s) {
    return { p.w / s, p.x / s, p.y / s, p.z / s };
}

inline Quat operator *(float s, Quat const& p) {
    return p * s;
}

inline Quat& operator +=(Quat& u, Quat const& s) { return u = u + s; }
inline Quat& operator -=(Quat& u, Quat const& s) { return u = u - s; }
inline Quat& operator *=(Quat& u, Quat const& s) { return u = u * s; }

inline Quat& operator *=(Quat& u, float s) { return u = u * s; }
inline Quat& operator /=(Quat& u, float s) { return u = u / s; }

//------------------------------------------------------------------------------
// Geometric functions
//

// Returns the rotation axis of the quaternion
// NOT normalized!
inline Vec3 RotationAxis(Quat const& q) {
    return { q.x, q.y, q.z };
}

// Returns the rotation angle of the quaternion
inline float RotationAngle(Quat const& q) {
    return 2.0f * Atan2(Length(q.vec()), q.w);
}

// Computes the dot-product of the two Quaternions
inline float Dot(Quat const& p, Quat const& q) {
    return p.w * q.w + p.x * q.x + p.y * q.y + p.z * q.z;
}

// Returns the squared norm of the quaternion
inline float LengthSquared(Quat const& q) {
    return Dot(q, q);
}

// Returns the length of the quaternion
inline float Length(Quat const& q) {
    return Sqrt(LengthSquared(q));
}

// Returns a normalized version of the given quaternion
inline Quat Normalize(Quat const& q) {
    return q / Length(q);
}

// Returns the conjugate of the Quat
inline Quat Conjugate(Quat const& q) { // reverse
    return { q.w, -q.x, -q.y, -q.z };
}

// Returns the inverse of the Quat
inline Quat Inverse(Quat const& q) {
    return Conjugate(q) / LengthSquared(q);
}

// Returns the angle between the two Quaternions
//
// NOTE:
// p and q must be normalized!
//
inline float Angle(Quat const& p, Quat const& q) {
    return 2.0f * Atan2(Length(p - q), Length(p + q));
}

// Performs a spherical linear interpolation
//
// NOTE:
// p and q must be normalized!
//
inline Quat Slerp(Quat const& p, Quat const& q, float t)
{
    const float phi = Angle(p, q);

    const float s = Sinc(phi);

    const float sp = (1.0f - t) * Sinc((1.0f - t) * phi) / s;
    const float sq = (       t) * Sinc((       t) * phi) / s;

    return Normalize(sp * p + sq * q);
}

// Performs a spherical linear interpolation
//
// NOTE:
// p and q must be normalized!
//
inline Quat SlerpShortest(Quat const& p, Quat const& q, float t) {
    return Slerp(p, Dot(p, q) >= 0.0f ? q : -q, t);
}

#if 0
// Quaternion integration step
// ^^  ^^         ^       ^
inline Quat Quergs(Quat const& dq /* = 1/2 q w dt */)
{
    // See:
    // http://www.cl.cam.ac.uk/techreports/UCAM-CL-TR-683.pdf (pp. 50)

    const auto s = Length(dq);

    return (Tan(s) / s) * dq;
}

inline Quat Quergs(Quat const& q, Quat const& omega, float dt) {
    return Normalize(q + Quergs((0.5f * dt) * (q * omega)));
}
#endif

//------------------------------------------------------------------------------
// Exponential functions
//

inline Quat Exp(Quat const& q)
{
    const auto v = Normalize(RotationAxis(q));
    const auto s = Length(v);

    return { Cos(s), Sin(s) * v };
}

inline Quat Log(Quat const& q)
{
    const auto v = Normalize(RotationAxis(q));
    const auto phi = RotationAngle(q);

    return { 0.0f, v / Sinc(phi) };
}

inline Quat Pow(Quat const& q, float t) {
    return Exp(t * Log(q));
}

#if 0
inline Quat Pow(Quat const& q, Quat const& p) {
    return Exp(p * Log(q));
}
#endif

//==============================================================================
// Quat
//==============================================================================

inline Euler::Euler(float heading, float elevation, float bank)
    : heading(heading)
    , elevation(elevation)
    , bank(bank)
{
}

inline Euler Euler::Identity() {
    return {0,0,0};
}

inline Euler Euler::Rotation(Mat3 const& m)
{
#if 0
    const float m20 = m(2,0);
    const float m01 = m(0,1);
    const float m11 = m(1,1);
    const float m21 = m(2,1);
    const float m22 = m(2,2);

    const float heading   = Atan2(-m20, m22);
    const float elevation = Asin(m21);
    const float bank      = Atan2(-m01, m11);
#endif
#if 0
    const float m20 = m(2,0);
    const float m01 = m(0,1);
    const float m11 = m(1,1);
    const float m21 = m(2,1);
    const float m22 = m(2,2);

    // Compute elevation using numerically more stable atan2 (instead of asin)

    const float heading   = Atan2(-m20, m22);
    const float elevation = Atan2(m21, Hypot(m01, m11));
    const float bank      = Atan2(-m01, m11);
#endif
#if 1
    const float m00 = m(0,0);
    const float m10 = m(1,0);
    const float m20 = m(2,0);
    const float m01 = m(0,1);
    const float m11 = m(1,1);
    const float m21 = m(2,1);
    const float m02 = m(0,2);
    const float m12 = m(1,2);
    const float m22 = m(2,2);

    const float heading   = Atan2(-m20, m22);
    const float elevation = Atan2(m21, Hypot(m01, m11));

    //
    // Solve for bank using M = Z * X * Y <==> Z = M * Y' * X', where X, Y and Z
    // are rotation matrices which rotate around the x-, y- and z-axis resp.
    //
    // Need to compute cos(heading) and sin(heading), where heading is
    // atan2(y,x) (see above).
    // Instead of directly computing cos and sin, use
    //
    //      cos(atan2(y,x)) = x/hypot(x,y)
    //      sin(atan2(y,x)) = y/hypot(x,y)
    //
    const float k = Hypot(m20, m22);
    const float c2 = -m20 / k;
    const float s2 =  m22 / k;

    const float bank = Atan2(m10 * c2 + m12 * s2, m00 * c2 + m02 * s2);
#endif

    return { heading, elevation, bank };
}

inline Euler Euler::Rotation(Quat const& q) {
    return Euler::Rotation(Mat3::Rotation(q));
}

} // namespace math
