/* surface_reconstructor.cpp
 * Reconstruct a triangle mesh from slice contours written by mesh_slicer.cpp.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

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

static double squared_distance(const Point& a, const Point& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
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

static std::vector<Slice> read_slices(const std::string& path)
{
    std::ifstream in(path.c_str());
    if( !in )
    {
        throw std::runtime_error("Invalid slice file.");
    }

    std::size_t nslices = 0;
    in >> nslices;
    if( !in )
    {
        throw std::runtime_error("Failed to read slice count.");
    }

    std::vector<Slice> slices(nslices);
    for( std::size_t t = 0; t != nslices; ++t )
    {
        std::size_t ncontours = 0;
        in >> ncontours;
        if( !in )
        {
            throw std::runtime_error("Failed to read contour count.");
        }

        slices[t].contours.reserve(ncontours);
        for( std::size_t c = 0; c != ncontours; ++c )
        {
            double area = 0.0;
            std::size_t npoints = 0;
            in >> area;
            in >> npoints;
            if( !in )
            {
                throw std::runtime_error("Failed to read contour header.");
            }

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

static std::vector<std::vector<Ring> > group_contours(std::vector<Slice>& slices)
{
    std::size_t max_contours = 0;
    for( std::size_t t = 0; t != slices.size(); ++t )
    {
        max_contours = std::max(max_contours, slices[t].contours.size());
    }

    std::vector<std::vector<Ring> > groups(max_contours);
    for( std::size_t contour = 0; contour != max_contours; ++contour )
    {
        for( std::size_t t = 0; t != slices.size(); ++t )
        {
            if( contour >= slices[t].contours.size() )
            {
                continue;
            }

            Ring ring = slices[t].contours[contour].ring;
            if( !groups[contour].empty() )
            {
                ring = best_aligned_ring(groups[contour].back(), ring);
            }
            groups[contour].push_back(ring);
        }
    }

    return groups;
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

static void add_cap(
    const std::vector<std::size_t>& ring_ids,
    const Ring& ring,
    bool reverse,
    std::vector<Point>& vertices,
    std::vector<Face>& faces)
{
    const std::size_t center_id = vertices.size();
    vertices.push_back(centroid(ring));

    for( std::size_t i = 0; i != ring_ids.size(); ++i )
    {
        const std::size_t next = (i + 1) % ring_ids.size();
        Face face;
        if( reverse )
        {
            face.vertices.push_back(center_id);
            face.vertices.push_back(ring_ids[next]);
            face.vertices.push_back(ring_ids[i]);
        }
        else
        {
            face.vertices.push_back(center_id);
            face.vertices.push_back(ring_ids[i]);
            face.vertices.push_back(ring_ids[next]);
        }
        faces.push_back(face);
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
        std::vector<std::vector<Ring> > contour_groups = group_contours(slices);

        std::vector<Point> vertices;
        std::vector<Face> faces;

        for( std::size_t contour = 0; contour != contour_groups.size(); ++contour )
        {
            const std::vector<Ring>& rings = contour_groups[contour];
            if( rings.size() < 2 )
            {
                continue;
            }

            std::vector<std::vector<std::size_t> > ring_ids;
            ring_ids.reserve(rings.size());
            for( std::size_t t = 0; t != rings.size(); ++t )
            {
                ring_ids.push_back(append_ring_vertices(rings[t], vertices));
            }

            add_cap(ring_ids.front(), rings.front(), true, vertices, faces);
            for( std::size_t t = 0; t + 1 != rings.size(); ++t )
            {
                add_surface_between(ring_ids[t], rings[t], ring_ids[t + 1], rings[t + 1], faces);
            }
            add_cap(ring_ids.back(), rings.back(), false, vertices, faces);
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
