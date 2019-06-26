/* mesh_slicer.cpp
 * Slice a mesh
 */

// https://www.cgal.org/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vector = K::Vector_3;
using Plane = K::Plane_3;
using Point = K::Point_3;
#include <CGAL/Surface_mesh.h>
using Mesh = CGAL::Surface_mesh<Point>;

// https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#PMPSlicer
#include <CGAL/Polygon_mesh_slicer.h>
using Slicer = CGAL::Polygon_mesh_slicer<Mesh, K>;

// https://doc.cgal.org/latest/Polygon/index.html
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
using Polyline = std::vector<Point>;

#include <CGAL/enum.h>

// standard C++ libraries
#include <iostream>
#include <fstream>
#include <cmath>

double calculate_area(Vector normal, Polyline poly3D, int i, std::vector<Polyline> other_lines)
{
    // calculate area: https://stackoverflow.com/a/30831214

    // ==== Find the "best" projection ====
    // i.e. the cartesian plane normal with lowest angle to the polygon normal

    Vector planes[3] = {Vector(1,0,0), Vector(0,1,0), Vector(0,0,1)};
    double absCosTheta = 1; // abs value of cos(theta), initially for 90deg
    int proj_dim;
    for( int dim = 0; dim != 3; ++dim )
    {
        Vector projection = planes[dim]; // already normalised
        Vector plane_normal = normal/sqrt(normal.squared_length());
        double absCosTheta_d = abs(projection*plane_normal); // [0,1]
        if( absCosTheta_d < absCosTheta ) // "closer" to this plane
        {
            // store
            absCosTheta = absCosTheta_d;
            proj_dim = dim;
        }
    }

    // ==== Now find the area ====

    double area;
    bool is_simple;
    bool is_negative = false;

    int removeDuplicateEnd = 1; //TODO: test if first and last are the same

    // compute 2D area of projection
    switch( proj_dim )
    {
        case 0: // [1,0,0]
            area = CGAL::polygon_area_2(poly3D.begin(),
                                        poly3D.end() -removeDuplicateEnd,
                                        CGAL::Projection_traits_yz_3<K>());
            is_simple = CGAL::is_simple_2(poly3D.begin(),
                                          poly3D.end() -removeDuplicateEnd,
                                          CGAL::Projection_traits_yz_3<K>());
            for( int j = 0; j != other_lines.size(); ++j )
            {
                if( j == i ) continue;
                Polyline other = other_lines.at(j);
                CGAL::Bounded_side side = CGAL::bounded_side_2(other.begin(),
                                                               other.end() -removeDuplicateEnd,
                                                               poly3D.at(0),
                                                               CGAL::Projection_traits_yz_3<K>());
                is_negative = ( side == CGAL::ON_BOUNDED_SIDE );
            }
            break;
        case 1: // [0,1,0]
            area = CGAL::polygon_area_2(poly3D.begin(),
                                        poly3D.end() -removeDuplicateEnd,
                                        CGAL::Projection_traits_xz_3<K>());
            is_simple = CGAL::is_simple_2(poly3D.begin(),
                                          poly3D.end() -removeDuplicateEnd,
                                          CGAL::Projection_traits_xz_3<K>());
            for( int j = 0; j != other_lines.size(); ++j )
            {
                if( j == i ) continue;
                Polyline other = other_lines.at(j);
                CGAL::Bounded_side side = CGAL::bounded_side_2(other.begin(),
                                                               other.end() -removeDuplicateEnd,
                                                               poly3D.at(0),
                                                               CGAL::Projection_traits_xz_3<K>());
                is_negative = ( side == CGAL::ON_BOUNDED_SIDE );
            }
            break;
        case 2: // [0,0,1]
            area = CGAL::polygon_area_2(poly3D.begin(),
                                        poly3D.end() -removeDuplicateEnd,
                                        CGAL::Projection_traits_xy_3<K>());
            is_simple = CGAL::is_simple_2(poly3D.begin(),
                                          poly3D.end() -removeDuplicateEnd,
                                          CGAL::Projection_traits_xy_3<K>());
            for( int j = 0; j != other_lines.size(); ++j )
            {
                if( j == i ) continue;
                Polyline other = other_lines.at(j);
                CGAL::Bounded_side side = CGAL::bounded_side_2(other.begin(),
                                                               other.end() -removeDuplicateEnd,
                                                               poly3D.at(0),
                                                               CGAL::Projection_traits_xy_3<K>());
                is_negative = ( side == CGAL::ON_BOUNDED_SIDE );
            }
            break;
    }

    if( !is_simple )
    {
        /* "simple"
         * = edges do not intersect, except consecutive edges at common vertex
         */
        std::cerr << "Warning! A polygon is not simple,"
                  << " its area may not be well-defined!" << std::endl;
    }

    // convert to 3D
    area = abs(area)/absCosTheta; // use abs! area can be negative if poly is cw
    if( is_negative ) area = -area;

    return area;
}

int main(int argc, char* argv[])
{
    /* Usage
     * $ mesh_slicer mesh_file.off centreline.txt outfile.txt
     * returns 0 on success, 1 on file error, 2 on mesh error
     */

    // construct slicer from mesh
    std::cout << "Reading mesh ..." << std::endl;
    std::ifstream meshstream(argv[1]); // supported formats: .off
    Mesh mesh;
    if( !meshstream )
    {
        std::cerr << "Invalid file!" << std::endl;
        return 1;
    }
    if( !(meshstream >> mesh) || mesh.is_empty() )
    {
        std::cerr << "Failed to read a mesh from file!" << std::endl;
        return 2;
    }
    if( !CGAL::is_triangle_mesh(mesh) )
    {
        std::cerr << "Only triangle meshes are supported!" << std::endl;
        return 2;
    }
    Slicer slicer(mesh);

    // read the centreline
    std::cout << "Reading centreline ..." << std::endl;
    std::ifstream linestream(argv[2]); // N rows of [Ox Oy Oz Nx Ny Nz]
    if( !linestream )
    {
        std::cerr << "Invalid centreline file!" << std::endl;
        return 1;
    }
    using Centreline = std::vector< std::tuple<Point, Vector> >;
    Centreline centreline;
    int ncoords;
    linestream >> ncoords;
    centreline.reserve(ncoords);
    for( int t = 0; t != ncoords; ++t )
    {
        double x, y, z; // dummy placeholders
        // origin
        linestream >> x;
        linestream >> y;
        linestream >> z;
        Point origin(x, y, z);
        // normal
        linestream >> x;
        linestream >> y;
        linestream >> z;
        Vector normal(x, y, z);
        // store
        centreline.emplace_back(std::make_tuple(origin, normal));
    }

    // slicing
    std::cout << "Slicing ..." << std::endl;
    std::ofstream outstream(argv[3]);
    if( !outstream )
    {
        std::cerr << "Invalid output file!" << std::endl;
        return 1;
    }
    outstream << centreline.size() << std::endl; // this many points
    for( std::size_t t = 0; t != centreline.size(); ++t )
    {
        Point  origin = std::get<0>( centreline.at(t) );
        Vector normal = std::get<1>( centreline.at(t) );

        // find polylines
        std::vector<Polyline> polylines; // there might be multiple
        slicer(Plane(origin, normal), std::back_inserter(polylines));
        outstream << polylines.size() << std::endl; // this many slices
        for( std::size_t i = 0; i != polylines.size(); ++i )
        {
            Polyline poly = polylines.at(i); // current

            // calculate area
            double area = calculate_area(normal, poly, i, polylines);
            outstream << area << std::endl;

            // write polygon
            outstream << poly.size() << std::endl;
            for( std::size_t j = 0; j != poly.size(); ++j )
            {
                outstream << poly.at(j) << std::endl;
            }
        }
        polylines.clear(); // needed?
    }

    std::cout << "Done!" << std::endl;
    return 0;
}
