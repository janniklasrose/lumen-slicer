/* tessellate_slices.cpp
 * Convert mesh_slicer slice output to a tessellated OFF mesh.
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2 = K::Point_2;

using VertexBase = CGAL::Triangulation_vertex_base_2<K>;
using FaceBase = CGAL::Delaunay_mesh_face_base_2<K>;
using TriangulationData = CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TriangulationData, CGAL::Exact_predicates_tag>;
using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;

struct Vec3
{
    double x;
    double y;
    double z;
};

struct Curve
{
    double area;
    std::vector<Vec3> points;
};

struct Slice
{
    std::vector<Curve> curves;
};

struct Frame
{
    Vec3 origin;
    Vec3 u;
    Vec3 v;
};

struct MeshOut
{
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t, 3> > faces;
};

static Vec3 operator+(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

static Vec3 operator-(const Vec3& a, const Vec3& b)
{
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

static Vec3 operator*(const Vec3& a, double s)
{
    return Vec3{a.x * s, a.y * s, a.z * s};
}

static double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static Vec3 cross(const Vec3& a, const Vec3& b)
{
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static double squared_length(const Vec3& a)
{
    return dot(a, a);
}

static Vec3 normalized(const Vec3& a)
{
    const double length = std::sqrt(squared_length(a));
    if(length == 0.0) throw std::runtime_error("zero-length vector");
    return a * (1.0 / length);
}

static bool same_point(const Vec3& a, const Vec3& b)
{
    const double tolerance = 1e-12;
    return squared_length(a - b) <= tolerance * tolerance;
}

static std::vector<Slice> read_slices(const std::string& path)
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

static std::vector<Vec3> unique_ring(const Curve& curve)
{
    std::vector<Vec3> ring = curve.points;
    while(ring.size() > 1 && same_point(ring.front(), ring.back()))
    {
        ring.pop_back();
    }
    return ring;
}

static Vec3 polygon_normal(const std::vector<Curve>& curves)
{
    Vec3 normal{0.0, 0.0, 0.0};
    for(const Curve& curve : curves)
    {
        const std::vector<Vec3> ring = unique_ring(curve);
        if(ring.size() < 3) continue;

        for(std::size_t i = 0; i != ring.size(); ++i)
        {
            const Vec3& current = ring[i];
            const Vec3& next = ring[(i + 1) % ring.size()];
            normal.x += (current.y - next.y) * (current.z + next.z);
            normal.y += (current.z - next.z) * (current.x + next.x);
            normal.z += (current.x - next.x) * (current.y + next.y);
        }
    }
    return normal;
}

static Frame make_frame(const Slice& slice)
{
    for(const Curve& curve : slice.curves)
    {
        const std::vector<Vec3> ring = unique_ring(curve);
        if(ring.size() < 3) continue;

        Vec3 normal = polygon_normal(slice.curves);
        if(squared_length(normal) == 0.0)
        {
            for(std::size_t i = 1; i + 1 < ring.size(); ++i)
            {
                normal = cross(ring[i] - ring[0], ring[i + 1] - ring[0]);
                if(squared_length(normal) > 0.0) break;
            }
        }
        normal = normalized(normal);

        Vec3 u{0.0, 0.0, 0.0};
        for(std::size_t i = 1; i != ring.size(); ++i)
        {
            u = ring[i] - ring[0];
            u = u - normal * dot(u, normal);
            if(squared_length(u) > 0.0) break;
        }
        u = normalized(u);
        const Vec3 v = normalized(cross(normal, u));

        return Frame{ring[0], u, v};
    }

    throw std::runtime_error("slice has no non-degenerate curve");
}

static Point2 project(const Frame& frame, const Vec3& p)
{
    const Vec3 offset = p - frame.origin;
    return Point2(dot(offset, frame.u), dot(offset, frame.v));
}

static Vec3 lift(const Frame& frame, const Point2& p)
{
    return frame.origin + frame.u * p.x() + frame.v * p.y();
}

static void insert_curve(CDT& cdt, const Frame& frame, const Curve& curve)
{
    const std::vector<Vec3> ring = unique_ring(curve);
    if(ring.size() < 3) return;

    std::vector<CDT::Vertex_handle> handles;
    handles.reserve(ring.size());
    for(const Vec3& point : ring)
    {
        handles.push_back(cdt.insert(project(frame, point)));
    }

    for(std::size_t i = 0; i != handles.size(); ++i)
    {
        cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
    }
}

static void append_domain_faces(const CDT& cdt, const Frame& frame, MeshOut& mesh)
{
    std::map<CDT::Vertex_handle, std::size_t> indices;

    for(CDT::Finite_faces_iterator face = cdt.finite_faces_begin();
        face != cdt.finite_faces_end();
        ++face)
    {
        if(!face->is_in_domain()) continue;

        std::array<std::size_t, 3> tri;
        for(int i = 0; i != 3; ++i)
        {
            const CDT::Vertex_handle vertex = face->vertex(i);
            std::map<CDT::Vertex_handle, std::size_t>::iterator found = indices.find(vertex);
            if(found == indices.end())
            {
                const std::size_t index = mesh.vertices.size();
                indices[vertex] = index;
                mesh.vertices.push_back(lift(frame, vertex->point()));
                tri[i] = index;
            }
            else
            {
                tri[i] = found->second;
            }
        }
        mesh.faces.push_back(tri);
    }
}

static void tessellate_slice(const Slice& slice, bool refine, MeshOut& mesh)
{
    if(slice.curves.empty()) return;

    const Frame frame = make_frame(slice);
    CDT cdt;
    for(const Curve& curve : slice.curves)
    {
        insert_curve(cdt, frame, curve);
    }

    if(cdt.dimension() != 2)
    {
        std::cerr << "Warning: skipping degenerate slice" << std::endl;
        return;
    }

    CGAL::mark_domain_in_triangulation(cdt);
    if(refine)
    {
        CGAL::refine_Delaunay_mesh_2(
            cdt,
            CGAL::parameters::domain_is_initialized(true).criteria(Criteria())
        );
    }

    append_domain_faces(cdt, frame, mesh);
}

static void write_off(const std::string& path, const MeshOut& mesh)
{
    std::ofstream out(path);
    if(!out) throw std::runtime_error("failed to open output file: " + path);

    out << "OFF\n";
    out << mesh.vertices.size() << " " << mesh.faces.size() << " 0\n";
    for(const Vec3& vertex : mesh.vertices)
    {
        out << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }
    for(const std::array<std::size_t, 3>& face : mesh.faces)
    {
        out << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
    }
}

static void usage(const char* program)
{
    std::cerr << "Usage: " << program << " slices.dat slices.off [--refine]\n";
}

int main(int argc, char* argv[])
{
    bool refine = false;
    std::vector<std::string> paths;

    for(int i = 1; i != argc; ++i)
    {
        const std::string arg(argv[i]);
        if(arg == "--refine")
        {
            refine = true;
        }
        else if(!arg.empty() && arg[0] == '-')
        {
            usage(argv[0]);
            std::cerr << "Unknown option: " << arg << std::endl;
            return 1;
        }
        else
        {
            paths.push_back(arg);
        }
    }

    if(paths.size() != 2)
    {
        usage(argv[0]);
        return 1;
    }

    try
    {
        const std::vector<Slice> slices = read_slices(paths[0]);
        MeshOut mesh;
        for(const Slice& slice : slices)
        {
            tessellate_slice(slice, refine, mesh);
        }
        write_off(paths[1], mesh);
    }
    catch(const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}
