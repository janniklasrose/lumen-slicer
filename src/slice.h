#ifndef LUMEN_SLICER_SLICE_H
#define LUMEN_SLICER_SLICE_H

#include "geometry.h"

#include <array>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct Curve
{
    double area;
    std::vector<Vec3> points;
};

struct Slice
{
    std::vector<Curve> curves;
};

struct MeshOut
{
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t, 3> > faces;
};

inline std::vector<Slice> read_slices(const std::string& path)
{
    std::ifstream in(path);
    if(!in) throw std::runtime_error("failed to open input file: " + path);

    int nslices = 0;
    if(!(in >> nslices) || nslices < 0) throw std::runtime_error("invalid slice count");

    std::vector<Slice> slices;
    slices.reserve(static_cast<std::size_t>(nslices));
    for(int i = 0; i != nslices; ++i)
    {
        int ncurves = 0;
        if(!(in >> ncurves) || ncurves < 0) throw std::runtime_error("invalid curve count");

        Slice slice;
        slice.curves.reserve(static_cast<std::size_t>(ncurves));
        for(int j = 0; j != ncurves; ++j)
        {
            Curve curve;
            int npoints = 0;
            if(!(in >> curve.area)) throw std::runtime_error("invalid curve area");
            if(!(in >> npoints) || npoints < 0) throw std::runtime_error("invalid point count");

            curve.points.reserve(static_cast<std::size_t>(npoints));
            for(int k = 0; k != npoints; ++k)
            {
                Vec3 p;
                if(!(in >> p.x >> p.y >> p.z)) throw std::runtime_error("invalid point coordinates");
                curve.points.push_back(p);
            }
            slice.curves.push_back(curve);
        }
        slices.push_back(slice);
    }

    return slices;
}

inline void write_off(const std::string& path, const MeshOut& mesh)
{
    std::ofstream out(path);
    if(!out) throw std::runtime_error("failed to open output file: " + path);

    out << "OFF\n";
    out << mesh.vertices.size() << " " << mesh.faces.size() << " 0\n";
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for(const Vec3& vertex : mesh.vertices)
    {
        out << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }
    for(const std::array<std::size_t, 3>& face : mesh.faces)
    {
        out << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
    }
}

#endif
