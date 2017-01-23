#pragma once

#define LIB_MATH_BOX_H 1

#include "Vector.h"
#include <cfloat>

namespace math {

//==============================================================================
// Box3
//==============================================================================

struct Box3
{
    Vec3 lower;
    Vec3 upper;

    Box3() = default;
    Box3(Vec3 const& lower, Vec3 const& upper);

    // Invalidates the box edges (sets lower > upper)
    void Invalidate();

    // Sets the bounding box to hold a single point
    void Clear(Vec3 const& point);

    // Returns whether the bounding box is invalid
    bool IsInvalid() const;

    // Returns whether the bounding box is valid
    bool IsValid() const;

    // Returns whether the bounding box is empty
    // Note: This includes invalid bounding boxes!
    bool IsEmpty() const;

    // Returns the center of the bounding box
    Vec3 Center() const;

    // Returns the size of the bounding box
    Vec3 Size() const;

    // Returns the size of the bounding box
    Vec3 SafeSize() const;

    // Insert the given point into this bounding box
    void Insert(Vec3 const& v);

    // Insert the given bounding box into this bounding box
    void Insert(Box3 const& box);
};

inline Box3::Box3(Vec3 const& lower, Vec3 const& upper)
    : lower(lower)
    , upper(upper)
{
}

inline void Box3::Invalidate()
{
    lower = Vec3( FLT_MAX);
    upper = Vec3(-FLT_MAX);
}

inline void Box3::Clear(Vec3 const& point)
{
    lower = point;
    upper = point;
}

inline bool Box3::IsInvalid() const {
    return upper.x < lower.x || upper.y < lower.y || upper.z < lower.z;
}

inline bool Box3::IsValid() const {
    return !IsInvalid();
}

inline bool Box3::IsEmpty() const {
    return upper.x <= lower.x || upper.y <= lower.y || upper.z <= lower.z;
}

inline Vec3 Box3::Center() const
{
//  return 0.5f * (lower + upper);
    return lower + 0.5f * (upper - lower);
}

inline Vec3 Box3::Size() const {
    return upper - lower;
}

inline Vec3 Box3::SafeSize() const
{
    auto s = Size();
    s.x = Max0(s.x);
    s.y = Max0(s.y);
    s.z = Max0(s.z);

    return s;
}

inline void Box3::Insert(Vec3 const& v)
{
    lower = Min(lower, v);
    upper = Max(upper, v);
}

inline void Box3::Insert(Box3 const& box)
{
    lower = Min(lower, box.lower);
    upper = Max(upper, box.upper);
}

//------------------------------------------------------------------------------
// Comparison operators
//

inline bool operator ==(Box3 const& lhs, Box3 const& rhs)
{
#if 0
    const bool lhs_invalid = lhs.IsInvalid();
    const bool rhs_invalid = rhs.IsInvalid();

    if (lhs_invalid == rhs_invalid)
        return true;

    if (lhs_invalid || rhs_invalid)
        return false;
#endif

    return lhs.lower == rhs.lower && lhs.upper == rhs.upper;
}

inline bool operator !=(Box3 const& lhs, Box3 const& rhs) {
    return !(lhs == rhs);
}

//------------------------------------------------------------------------------
// Arithmetic
//

inline Box3 operator +(Box3 const& box, Vec3 const& d) { return { box.lower + d, box.upper + d }; }
inline Box3 operator -(Box3 const& box, Vec3 const& d) { return { box.lower - d, box.upper - d }; }
inline Box3 operator *(Box3 const& box, Vec3 const& d) { return { box.lower * d, box.upper * d }; }
inline Box3 operator /(Box3 const& box, Vec3 const& d) { return { box.lower / d, box.upper / d }; }

inline Box3 operator +(Vec3 const& d, Box3 const& box) { return { d + box.lower, d + box.upper }; }
inline Box3 operator -(Vec3 const& d, Box3 const& box) { return { d - box.lower, d - box.upper }; }
inline Box3 operator *(Vec3 const& d, Box3 const& box) { return { d * box.lower, d * box.upper }; }

inline Box3& operator +=(Box3& r, Vec3 const& d) { r = r + d; return r; }
inline Box3& operator -=(Box3& r, Vec3 const& d) { r = r - d; return r; }
inline Box3& operator *=(Box3& r, Vec3 const& s) { r = r * s; return r; }
inline Box3& operator /=(Box3& r, Vec3 const& s) { r = r / s; return r; }

//------------------------------------------------------------------------------
// Geometric functions
//

// Computes the intersection of two bounding boxes
inline Box3 Intersect(Box3 const& lhs, Box3 const& rhs) {
    return { Max(lhs.lower, rhs.lower), Min(lhs.upper, rhs.upper) };
}

// Computes the union of two bounding boxes
inline Box3 Merge(Box3 const& lhs, Box3 const& rhs) {
    return { Min(lhs.lower, rhs.lower), Max(lhs.upper, rhs.upper) };
}

// Computes the union of a bounding box and a point
inline Box3 Merge(Box3 const& lhs, Vec3 const& rhs) {
    return { Min(lhs.lower, rhs), Max(lhs.upper, rhs) };
}

// Computes the union of a bounding box and a point
inline Box3 Merge(Vec3 const& lhs, Box3 const& rhs) { return Merge(rhs, lhs); }

// Returns the volume of the bounding box
inline float Volume(Box3 const& box) { return ReduceMul(box.Size()); }

// Returns the surface area, i.e. the volume of the boundary.
inline float SurfaceArea(Box3 const& box)
{
    const auto s = box.Size();
    const auto t = s.x * s.y + s.y * s.z + s.z * s.x;

    return 2.0f * t;
}

// Returns the surface area, i.e. the volume of the boundary.
inline float SafeSurfaceArea(Box3 const& box)
{
    const auto s = box.SafeSize();
    const auto t = s.x * s.y + s.y * s.z + s.z * s.x;

    return 2.0f * t;
}

// Returns the vertices of the box
inline void ComputeVertices(Box3 const& box, Vec3 vertices[8])
{
    //
    // 0 = min
    // 6 = max
    //                 3 ---- 2
    //     y          /|    / |
    //     |        7 ---- 6  |
    //     +-- x    |  0 --|- 1
    //    /         | /    | /
    //   z          4 ---- 5
    //

    const float minx = box.lower.x;
    const float miny = box.lower.y;
    const float minz = box.lower.z;
    const float maxx = box.upper.x;
    const float maxy = box.upper.y;
    const float maxz = box.upper.z;

    vertices[0] = Vec3(minx, miny, minz);
    vertices[1] = Vec3(maxx, miny, minz);
    vertices[2] = Vec3(maxx, maxy, minz);
    vertices[3] = Vec3(minx, maxy, minz);
    vertices[4] = Vec3(minx, miny, maxz);
    vertices[5] = Vec3(maxx, miny, maxz);
    vertices[6] = Vec3(maxx, maxy, maxz);
    vertices[7] = Vec3(minx, maxy, maxz);
}

} // namespace math
