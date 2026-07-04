/* surface_reconstructor.cpp
 * Reconstruct a triangle mesh from slice contours written by mesh_slicer.cpp.
 */

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2 = K::Point_2;
using VertexBase = CGAL::Triangulation_vertex_base_2<K>;
using FaceBase = CGAL::Delaunay_mesh_face_base_2<K>;
using TriangulationData = CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TriangulationData, CGAL::Exact_predicates_tag>;

struct Point
{
    double x;
    double y;
    double z;
};

struct Face
{
    std::vector<std::size_t> vertices;
};

using Ring = std::vector<Point>;

struct Contour
{
    double area;
    Ring ring;
};

struct Slice
{
    std::vector<Contour> contours;
};

struct ContourSample
{
    std::size_t slice_index;
    double area;
    Ring ring;
};

using ContourTrack = std::vector<ContourSample>;

struct CapContour
{
    std::vector<std::size_t> ring_ids;
    Ring ring;
};

struct Frame
{
    Point origin;
    Point u;
    Point v;
};

static Point operator+(const Point& a, const Point& b)
{
    Point result = {a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

static Point operator-(const Point& a, const Point& b)
{
    Point result = {a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

static Point operator*(const Point& a, double scale)
{
    Point result = {a.x * scale, a.y * scale, a.z * scale};
    return result;
}

static double dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static Point cross(const Point& a, const Point& b)
{
    Point result = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return result;
}

static double squared_distance(const Point& a, const Point& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

static double squared_length(const Point& a)
{
    return dot(a, a);
}

static Point normalized(const Point& a)
{
    const double length = std::sqrt(squared_length(a));
    if( length == 0.0 )
    {
        throw std::runtime_error("Zero-length vector.");
    }
    return a * (1.0 / length);
}

static bool same_point(const Point& a, const Point& b)
{
    const double eps = 1e-12;
    return squared_distance(a, b) <= eps * eps;
}

static void remove_duplicate_endpoint(Ring& ring)
{
    if( ring.size() > 1 && same_point(ring.front(), ring.back()) )
    {
        ring.pop_back();
    }
}

static std::size_t read_count(std::istream& in, const std::string& description)
{
    long long value = 0;
    if( !(in >> value) || value < 0 )
    {
        throw std::runtime_error("Invalid " + description + ".");
    }
    return static_cast<std::size_t>(value);
}

static Point centroid(const Ring& ring)
{
    Point c = {0.0, 0.0, 0.0};
    for( std::size_t i = 0; i != ring.size(); ++i )
    {
        c.x += ring[i].x;
        c.y += ring[i].y;
        c.z += ring[i].z;
    }
    const double n = static_cast<double>(ring.size());
    c.x /= n;
    c.y /= n;
    c.z /= n;
    return c;
}

static Point ring_normal(const Ring& ring)
{
    Point normal = {0.0, 0.0, 0.0};
    for( std::size_t i = 0; i != ring.size(); ++i )
    {
        const Point& current = ring[i];
        const Point& next = ring[(i + 1) % ring.size()];
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }
    return normal;
}

static Point non_collinear_normal(const Ring& ring)
{
    Point normal = ring_normal(ring);
    if( squared_length(normal) != 0.0 )
    {
        return normal;
    }

    for( std::size_t i = 1; i + 1 < ring.size(); ++i )
    {
        normal = cross(ring[i] - ring[0], ring[i + 1] - ring[0]);
        if( squared_length(normal) != 0.0 )
        {
            return normal;
        }
    }
    Point zero = {0.0, 0.0, 0.0};
    return zero;
}

static bool make_frame(const std::vector<CapContour>& contours, Frame& frame)
{
    for( std::size_t contour = 0; contour != contours.size(); ++contour )
    {
        const Ring& ring = contours[contour].ring;
        if( ring.size() < 3 )
        {
            continue;
        }

        Point normal = non_collinear_normal(ring);
        if( squared_length(normal) == 0.0 )
        {
            continue;
        }
        normal = normalized(normal);

        Point u = {0.0, 0.0, 0.0};
        for( std::size_t i = 1; i != ring.size(); ++i )
        {
            u = ring[i] - ring[0];
            u = u - normal * dot(u, normal);
            if( squared_length(u) > 0.0 )
            {
                break;
            }
        }
        if( squared_length(u) == 0.0 )
        {
            continue;
        }
        u = normalized(u);

        frame.origin = ring[0];
        frame.u = u;
        frame.v = normalized(cross(normal, u));
        return true;
    }

    return false;
}

static Point2 project(const Frame& frame, const Point& point)
{
    const Point offset = point - frame.origin;
    return Point2(dot(offset, frame.u), dot(offset, frame.v));
}

static Point lift(const Frame& frame, const Point2& point)
{
    return frame.origin + frame.u * point.x() + frame.v * point.y();
}

static Ring shifted_ring(const Ring& ring, std::size_t offset, bool reversed)
{
    Ring shifted;
    shifted.reserve(ring.size());

    if( reversed )
    {
        for( std::size_t k = 0; k != ring.size(); ++k )
        {
            const std::size_t idx = (offset + ring.size() - k) % ring.size();
            shifted.push_back(ring[idx]);
        }
    }
    else
    {
        for( std::size_t k = 0; k != ring.size(); ++k )
        {
            const std::size_t idx = (offset + k) % ring.size();
            shifted.push_back(ring[idx]);
        }
    }

    return shifted;
}

static double alignment_score(const Ring& a, const Ring& b, std::size_t offset, bool reversed)
{
    const std::size_t comparisons = std::min(a.size(), b.size());
    double score = 0.0;

    for( std::size_t k = 0; k != comparisons; ++k )
    {
        const std::size_t bidx = reversed
            ? (offset + b.size() - k) % b.size()
            : (offset + k) % b.size();
        const std::size_t aidx = (k * a.size()) / comparisons;
        score += squared_distance(a[aidx], b[bidx]);
    }

    return score;
}

static Ring best_aligned_ring(const Ring& previous, const Ring& current)
{
    double best_score = std::numeric_limits<double>::max();
    std::size_t best_offset = 0;
    bool best_reversed = false;

    for( std::size_t offset = 0; offset != current.size(); ++offset )
    {
        for( int reverse = 0; reverse != 2; ++reverse )
        {
            const bool reversed = ( reverse == 1 );
            const double score = alignment_score(previous, current, offset, reversed);
            if( score < best_score )
            {
                best_score = score;
                best_offset = offset;
                best_reversed = reversed;
            }
        }
    }

    return shifted_ring(current, best_offset, best_reversed);
}

static int area_class(double area)
{
    return area < 0.0 ? -1 : 1;
}

static double contour_match_score(const ContourSample& previous, const Contour& current)
{
    const Point previous_center = centroid(previous.ring);
    const Point current_center = centroid(current.ring);
    const double previous_radius = std::sqrt(std::fabs(previous.area));
    const double current_radius = std::sqrt(std::fabs(current.area));
    const double scale = std::max(std::max(previous_radius, current_radius), 1e-12);
    const double center_score = squared_distance(previous_center, current_center) / (scale * scale);
    const double radius_delta = previous_radius - current_radius;
    const double area_score = (radius_delta * radius_delta) / (scale * scale);
    return center_score + area_score;
}

static std::vector<Slice> read_slices(const std::string& path)
{
    std::ifstream in(path.c_str());
    if( !in )
    {
        throw std::runtime_error("Invalid slice file.");
    }

    const std::size_t nslices = read_count(in, "slice count");

    std::vector<Slice> slices(nslices);
    for( std::size_t t = 0; t != nslices; ++t )
    {
        const std::size_t ncontours = read_count(in, "contour count");

        slices[t].contours.reserve(ncontours);
        for( std::size_t c = 0; c != ncontours; ++c )
        {
            double area = 0.0;
            in >> area;
            if( !in )
            {
                throw std::runtime_error("Invalid contour area.");
            }
            const std::size_t npoints = read_count(in, "point count");

            Contour contour;
            contour.area = area;
            contour.ring.reserve(npoints);
            for( std::size_t i = 0; i != npoints; ++i )
            {
                Point p;
                in >> p.x >> p.y >> p.z;
                if( !in )
                {
                    throw std::runtime_error("Failed to read contour point.");
                }
                contour.ring.push_back(p);
            }

            remove_duplicate_endpoint(contour.ring);
            if( contour.ring.size() < 3 )
            {
                std::ostringstream msg;
                msg << "Contour " << c << " in slice " << t << " has fewer than 3 points.";
                throw std::runtime_error(msg.str());
            }

            slices[t].contours.push_back(contour);
        }
    }

    return slices;
}

static std::vector<ContourTrack> group_contours(const std::vector<Slice>& slices)
{
    struct Candidate
    {
        double score;
        std::size_t active_index;
        std::size_t contour_index;
    };

    std::vector<ContourTrack> tracks;
    std::vector<std::size_t> active_tracks;

    for( std::size_t t = 0; t != slices.size(); ++t )
    {
        const std::vector<Contour>& contours = slices[t].contours;
        std::vector<Candidate> candidates;
        for( std::size_t active = 0; active != active_tracks.size(); ++active )
        {
            const ContourSample& previous = tracks[active_tracks[active]].back();
            for( std::size_t contour = 0; contour != contours.size(); ++contour )
            {
                if( area_class(previous.area) != area_class(contours[contour].area) )
                {
                    continue;
                }
                Candidate candidate;
                candidate.score = contour_match_score(previous, contours[contour]);
                candidate.active_index = active;
                candidate.contour_index = contour;
                candidates.push_back(candidate);
            }
        }

        std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b)
            {
                return a.score < b.score;
            });

        std::vector<bool> used_active(active_tracks.size(), false);
        std::vector<bool> used_contours(contours.size(), false);
        std::vector<std::size_t> next_active_tracks;

        for( std::size_t i = 0; i != candidates.size(); ++i )
        {
            const Candidate& candidate = candidates[i];
            if( used_active[candidate.active_index] || used_contours[candidate.contour_index] )
            {
                continue;
            }

            const std::size_t track_index = active_tracks[candidate.active_index];
            const ContourSample& previous = tracks[track_index].back();
            const Contour& contour = contours[candidate.contour_index];

            ContourSample sample;
            sample.slice_index = t;
            sample.area = contour.area;
            sample.ring = best_aligned_ring(previous.ring, contour.ring);
            tracks[track_index].push_back(sample);

            used_active[candidate.active_index] = true;
            used_contours[candidate.contour_index] = true;
            next_active_tracks.push_back(track_index);
        }

        for( std::size_t contour = 0; contour != contours.size(); ++contour )
        {
            if( used_contours[contour] )
            {
                continue;
            }

            ContourTrack track;
            ContourSample sample;
            sample.slice_index = t;
            sample.area = contours[contour].area;
            sample.ring = contours[contour].ring;
            track.push_back(sample);
            tracks.push_back(track);
            next_active_tracks.push_back(tracks.size() - 1);
        }

        active_tracks = next_active_tracks;
    }

    return tracks;
}

static std::vector<std::size_t> append_ring_vertices(
    const Ring& ring,
    std::vector<Point>& vertices)
{
    std::vector<std::size_t> ids;
    ids.reserve(ring.size());
    for( std::size_t i = 0; i != ring.size(); ++i )
    {
        ids.push_back(vertices.size());
        vertices.push_back(ring[i]);
    }
    return ids;
}

static void insert_cap_contour(
    CDT& cdt,
    const Frame& frame,
    const CapContour& contour,
    std::map<CDT::Vertex_handle, std::size_t>& vertex_ids)
{
    std::vector<CDT::Vertex_handle> handles;
    handles.reserve(contour.ring.size());

    for( std::size_t i = 0; i != contour.ring.size(); ++i )
    {
        const CDT::Vertex_handle handle = cdt.insert(project(frame, contour.ring[i]));
        handles.push_back(handle);
        vertex_ids[handle] = contour.ring_ids[i];
    }

    for( std::size_t i = 0; i != handles.size(); ++i )
    {
        cdt.insert_constraint(handles[i], handles[(i + 1) % handles.size()]);
    }
}

static std::size_t cap_vertex_id(
    const CDT::Vertex_handle& vertex,
    const Frame& frame,
    std::map<CDT::Vertex_handle, std::size_t>& vertex_ids,
    std::vector<Point>& vertices)
{
    std::map<CDT::Vertex_handle, std::size_t>::const_iterator found = vertex_ids.find(vertex);
    if( found != vertex_ids.end() )
    {
        return found->second;
    }

    const std::size_t vertex_id = vertices.size();
    vertices.push_back(lift(frame, vertex->point()));
    vertex_ids[vertex] = vertex_id;
    return vertex_id;
}

static void add_cap(
    const std::vector<CapContour>& contours,
    bool reverse,
    std::vector<Point>& vertices,
    std::vector<Face>& faces)
{
    if( contours.empty() )
    {
        return;
    }

    Frame frame;
    if( !make_frame(contours, frame) )
    {
        throw std::runtime_error("Failed to construct cap frame.");
    }

    CDT cdt;
    std::map<CDT::Vertex_handle, std::size_t> vertex_ids;
    for( std::size_t contour = 0; contour != contours.size(); ++contour )
    {
        insert_cap_contour(cdt, frame, contours[contour], vertex_ids);
    }

    if( cdt.dimension() != 2 )
    {
        throw std::runtime_error("Failed to triangulate cap.");
    }

    CGAL::mark_domain_in_triangulation(cdt);
    for( CDT::Finite_faces_iterator face = cdt.finite_faces_begin();
         face != cdt.finite_faces_end();
         ++face )
    {
        if( !face->is_in_domain() )
        {
            continue;
        }

        Face output;
        if( reverse )
        {
            output.vertices.push_back(cap_vertex_id(face->vertex(0), frame, vertex_ids, vertices));
            output.vertices.push_back(cap_vertex_id(face->vertex(2), frame, vertex_ids, vertices));
            output.vertices.push_back(cap_vertex_id(face->vertex(1), frame, vertex_ids, vertices));
        }
        else
        {
            output.vertices.push_back(cap_vertex_id(face->vertex(0), frame, vertex_ids, vertices));
            output.vertices.push_back(cap_vertex_id(face->vertex(1), frame, vertex_ids, vertices));
            output.vertices.push_back(cap_vertex_id(face->vertex(2), frame, vertex_ids, vertices));
        }
        faces.push_back(output);
    }
}

static void add_surface_between(
    const std::vector<std::size_t>& a_ids,
    const Ring& a,
    const std::vector<std::size_t>& b_ids,
    const Ring& b,
    std::vector<Face>& faces)
{
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t advanced_a = 0;
    std::size_t advanced_b = 0;

    while( advanced_a != a.size() || advanced_b != b.size() )
    {
        const std::size_t inext = (i + 1) % a.size();
        const std::size_t jnext = (j + 1) % b.size();
        const bool must_advance_a = ( advanced_b == b.size() );
        const bool must_advance_b = ( advanced_a == a.size() );

        if( must_advance_a ||
            ( !must_advance_b && squared_distance(a[inext], b[j]) <= squared_distance(a[i], b[jnext]) ) )
        {
            Face face;
            face.vertices.push_back(a_ids[i]);
            face.vertices.push_back(a_ids[inext]);
            face.vertices.push_back(b_ids[j]);
            faces.push_back(face);
            i = inext;
            ++advanced_a;
        }
        else
        {
            Face face;
            face.vertices.push_back(a_ids[i]);
            face.vertices.push_back(b_ids[j]);
            face.vertices.push_back(b_ids[jnext]);
            faces.push_back(face);
            j = jnext;
            ++advanced_b;
        }
    }
}

static void write_off(
    const std::string& path,
    const std::vector<Point>& vertices,
    const std::vector<Face>& faces)
{
    std::ofstream out(path.c_str());
    if( !out )
    {
        throw std::runtime_error("Invalid output file.");
    }

    out << "OFF\n";
    out << vertices.size() << " " << faces.size() << " 0\n";
    for( std::size_t i = 0; i != vertices.size(); ++i )
    {
        out << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << "\n";
    }
    for( std::size_t i = 0; i != faces.size(); ++i )
    {
        out << faces[i].vertices.size();
        for( std::size_t j = 0; j != faces[i].vertices.size(); ++j )
        {
            out << " " << faces[i].vertices[j];
        }
        out << "\n";
    }
}

int main(int argc, char* argv[])
{
    /* Usage
     * $ surface_reconstructor slices.dat mesh.off
     * returns 0 on success, 1 on file/format error, 2 on reconstruction error
     */
    if( argc != 3 )
    {
        std::cerr << "Usage: " << argv[0] << " slices.dat mesh.off" << std::endl;
        return 1;
    }

    try
    {
        std::vector<Slice> slices = read_slices(argv[1]);
        std::vector<ContourTrack> contour_tracks = group_contours(slices);

        std::vector<Point> vertices;
        std::vector<Face> faces;
        std::map<std::size_t, std::vector<CapContour> > start_caps;
        std::map<std::size_t, std::vector<CapContour> > end_caps;

        for( std::size_t contour = 0; contour != contour_tracks.size(); ++contour )
        {
            const ContourTrack& track = contour_tracks[contour];
            if( track.size() < 2 )
            {
                continue;
            }

            std::vector<std::vector<std::size_t> > ring_ids;
            ring_ids.reserve(track.size());
            for( std::size_t t = 0; t != track.size(); ++t )
            {
                ring_ids.push_back(append_ring_vertices(track[t].ring, vertices));
            }

            CapContour front_cap;
            front_cap.ring_ids = ring_ids.front();
            front_cap.ring = track.front().ring;
            start_caps[track.front().slice_index].push_back(front_cap);
            for( std::size_t t = 0; t + 1 != track.size(); ++t )
            {
                add_surface_between(ring_ids[t], track[t].ring, ring_ids[t + 1], track[t + 1].ring, faces);
            }
            CapContour back_cap;
            back_cap.ring_ids = ring_ids.back();
            back_cap.ring = track.back().ring;
            end_caps[track.back().slice_index].push_back(back_cap);
        }

        for( std::map<std::size_t, std::vector<CapContour> >::const_iterator cap = start_caps.begin();
             cap != start_caps.end();
             ++cap )
        {
            add_cap(cap->second, true, vertices, faces);
        }
        for( std::map<std::size_t, std::vector<CapContour> >::const_iterator cap = end_caps.begin();
             cap != end_caps.end();
             ++cap )
        {
            add_cap(cap->second, false, vertices, faces);
        }

        if( faces.empty() )
        {
            throw std::runtime_error("At least two slices with matching contour indices are required.");
        }

        write_off(argv[2], vertices, faces);
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 2;
    }

    return 0;
}
