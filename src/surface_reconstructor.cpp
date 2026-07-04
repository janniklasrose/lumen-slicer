/* surface_reconstructor.cpp
 * Reconstruct a triangle mesh from slice contours written by mesh_slicer.cpp.
 */

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include "geometry.h"
#include "slice.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
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

using Ring = std::vector<Vec3>;

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
    Vec3 origin;
    Vec3 u;
    Vec3 v;
};

struct FormatError : public std::runtime_error
{
    explicit FormatError(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

static double squared_distance(const Vec3& a, const Vec3& b)
{
    return squared_length(a - b);
}

static Vec3 centroid(const Ring& ring)
{
    Vec3 c = {0.0, 0.0, 0.0};
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

static Ring unique_ring(const Curve& curve)
{
    Ring ring = curve.points;
    while( ring.size() > 1 && same_point(ring.front(), ring.back(), 1e-12) )
    {
        ring.pop_back();
    }
    return ring;
}

static Vec3 ring_normal(const Ring& ring)
{
    Vec3 normal = {0.0, 0.0, 0.0};
    for( std::size_t i = 0; i != ring.size(); ++i )
    {
        const Vec3& current = ring[i];
        const Vec3& next = ring[(i + 1) % ring.size()];
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }
    return normal;
}

static Vec3 non_collinear_normal(const Ring& ring)
{
    Vec3 normal = ring_normal(ring);
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
    Vec3 zero = {0.0, 0.0, 0.0};
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

        Vec3 normal = non_collinear_normal(ring);
        if( squared_length(normal) == 0.0 )
        {
            continue;
        }
        normal = normalized(normal);

        Vec3 u = {0.0, 0.0, 0.0};
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

static Point2 project(const Frame& frame, const Vec3& point)
{
    const Vec3 offset = point - frame.origin;
    return Point2(dot(offset, frame.u), dot(offset, frame.v));
}

static Vec3 lift(const Frame& frame, const Point2& point)
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

static ContourSample make_contour_sample(
    std::size_t slice_index,
    std::size_t contour_index,
    const Curve& curve)
{
    ContourSample sample;
    sample.slice_index = slice_index;
    sample.area = curve.area;
    sample.ring = unique_ring(curve);

    if( sample.ring.size() < 3 )
    {
        std::ostringstream msg;
        msg << "Contour " << contour_index << " in slice " << slice_index
            << " has fewer than 3 points.";
        throw FormatError(msg.str());
    }

    return sample;
}

static double contour_match_score(const ContourSample& previous, const ContourSample& current)
{
    const Vec3 previous_center = centroid(previous.ring);
    const Vec3 current_center = centroid(current.ring);
    const double previous_radius = std::sqrt(std::fabs(previous.area));
    const double current_radius = std::sqrt(std::fabs(current.area));
    const double scale = std::max(std::max(previous_radius, current_radius), 1e-12);
    const double center_score = squared_distance(previous_center, current_center) / (scale * scale);
    const double radius_delta = previous_radius - current_radius;
    const double area_score = (radius_delta * radius_delta) / (scale * scale);
    return center_score + area_score;
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
        const std::vector<Curve>& curves = slices[t].curves;
        std::vector<ContourSample> samples;
        samples.reserve(curves.size());
        for( std::size_t contour = 0; contour != curves.size(); ++contour )
        {
            samples.push_back(make_contour_sample(t, contour, curves[contour]));
        }

        std::vector<Candidate> candidates;
        for( std::size_t active = 0; active != active_tracks.size(); ++active )
        {
            const ContourSample& previous = tracks[active_tracks[active]].back();
            for( std::size_t contour = 0; contour != samples.size(); ++contour )
            {
                if( area_class(previous.area) != area_class(samples[contour].area) )
                {
                    continue;
                }
                Candidate candidate;
                candidate.score = contour_match_score(previous, samples[contour]);
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
        std::vector<bool> used_contours(curves.size(), false);
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
            const ContourSample& candidate_sample = samples[candidate.contour_index];

            ContourSample sample;
            sample.slice_index = t;
            sample.area = candidate_sample.area;
            sample.ring = best_aligned_ring(previous.ring, candidate_sample.ring);
            tracks[track_index].push_back(sample);

            used_active[candidate.active_index] = true;
            used_contours[candidate.contour_index] = true;
            next_active_tracks.push_back(track_index);
        }

        for( std::size_t contour = 0; contour != samples.size(); ++contour )
        {
            if( used_contours[contour] )
            {
                continue;
            }

            ContourTrack track;
            track.push_back(samples[contour]);
            tracks.push_back(track);
            next_active_tracks.push_back(tracks.size() - 1);
        }

        active_tracks = next_active_tracks;
    }

    return tracks;
}

static std::vector<std::size_t> append_ring_vertices(
    const Ring& ring,
    MeshOut& mesh)
{
    std::vector<std::size_t> ids;
    ids.reserve(ring.size());
    for( std::size_t i = 0; i != ring.size(); ++i )
    {
        ids.push_back(mesh.vertices.size());
        mesh.vertices.push_back(ring[i]);
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
    MeshOut& mesh)
{
    std::map<CDT::Vertex_handle, std::size_t>::const_iterator found = vertex_ids.find(vertex);
    if( found != vertex_ids.end() )
    {
        return found->second;
    }

    const std::size_t vertex_id = mesh.vertices.size();
    mesh.vertices.push_back(lift(frame, vertex->point()));
    vertex_ids[vertex] = vertex_id;
    return vertex_id;
}

static void add_cap(
    const std::vector<CapContour>& contours,
    bool reverse,
    MeshOut& mesh)
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

        std::array<std::size_t, 3> output;
        if( reverse )
        {
            output = {{
                cap_vertex_id(face->vertex(0), frame, vertex_ids, mesh),
                cap_vertex_id(face->vertex(2), frame, vertex_ids, mesh),
                cap_vertex_id(face->vertex(1), frame, vertex_ids, mesh)
            }};
        }
        else
        {
            output = {{
                cap_vertex_id(face->vertex(0), frame, vertex_ids, mesh),
                cap_vertex_id(face->vertex(1), frame, vertex_ids, mesh),
                cap_vertex_id(face->vertex(2), frame, vertex_ids, mesh)
            }};
        }
        mesh.faces.push_back(output);
    }
}

static void add_surface_between(
    const std::vector<std::size_t>& a_ids,
    const Ring& a,
    const std::vector<std::size_t>& b_ids,
    const Ring& b,
    MeshOut& mesh)
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
            mesh.faces.push_back({{a_ids[i], a_ids[inext], b_ids[j]}});
            i = inext;
            ++advanced_a;
        }
        else
        {
            mesh.faces.push_back({{a_ids[i], b_ids[j], b_ids[jnext]}});
            j = jnext;
            ++advanced_b;
        }
    }
}

static std::vector<Slice> load_slices(const std::string& path)
{
    try
    {
        return read_slices(path);
    }
    catch( const std::exception& e )
    {
        throw FormatError(e.what());
    }
}

static void write_mesh(const std::string& path, const MeshOut& mesh)
{
    try
    {
        write_off(path, mesh);
    }
    catch( const std::exception& e )
    {
        throw FormatError(e.what());
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
        std::vector<Slice> slices = load_slices(argv[1]);
        std::vector<ContourTrack> contour_tracks = group_contours(slices);

        MeshOut mesh;
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
                ring_ids.push_back(append_ring_vertices(track[t].ring, mesh));
            }

            CapContour front_cap;
            front_cap.ring_ids = ring_ids.front();
            front_cap.ring = track.front().ring;
            start_caps[track.front().slice_index].push_back(front_cap);
            for( std::size_t t = 0; t + 1 != track.size(); ++t )
            {
                add_surface_between(ring_ids[t], track[t].ring, ring_ids[t + 1], track[t + 1].ring, mesh);
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
            add_cap(cap->second, true, mesh);
        }
        for( std::map<std::size_t, std::vector<CapContour> >::const_iterator cap = end_caps.begin();
             cap != end_caps.end();
             ++cap )
        {
            add_cap(cap->second, false, mesh);
        }

        if( mesh.faces.empty() )
        {
            throw std::runtime_error("At least two slices with matching contour indices are required.");
        }

        write_mesh(argv[2], mesh);
    }
    catch( const FormatError& e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 2;
    }

    return 0;
}
