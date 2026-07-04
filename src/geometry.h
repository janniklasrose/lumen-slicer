#ifndef LUMEN_SLICER_GEOMETRY_H
#define LUMEN_SLICER_GEOMETRY_H

#include <cmath>
#include <stdexcept>

struct Vec3
{
    double x;
    double y;
    double z;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3 operator*(const Vec3& a, double s)
{
    return Vec3{a.x * s, a.y * s, a.z * s};
}

inline double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

inline double squared_length(const Vec3& a)
{
    return dot(a, a);
}

inline Vec3 normalized(const Vec3& a)
{
    const double length = std::sqrt(squared_length(a));
    if(length == 0.0) throw std::runtime_error("zero-length vector");
    return a * (1.0 / length);
}

inline bool same_point(const Vec3& a, const Vec3& b, double tolerance)
{
    return squared_length(a - b) <= tolerance * tolerance;
}

#endif
