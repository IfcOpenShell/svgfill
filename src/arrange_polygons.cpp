// #define SVGFILL_DEBUG
// #define SVGFILL_MAIN

#ifndef SVGFILL_MAIN
#include "svgfill.h"
#endif

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/create_offset_polygons_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/box_intersection_d.h>

#include <vector>
#include <iostream>

#include "graph_2d.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef std::vector<Polygon_2> Polygon_list;
typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef typename Traits_2::Point_2 Arr_Point_2;
typedef typename Traits_2::X_monotone_curve_2 Arr_Segment_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;

template <typename K>
using Triangle = std::array<CGAL::Point_2<K>, 3>;

template <typename K2>
Polygon_2 convert_polygon(const CGAL::Polygon_2<K2>& poly) {
    Polygon_2 exact_poly;
    typedef CGAL::Cartesian_converter<K2, K> Converter_Epick_to_Epeck;
    Converter_Epick_to_Epeck converter;
    for (auto vit = poly.vertices_begin(); vit != poly.vertices_end(); ++vit) {
        exact_poly.push_back(converter(*vit));  // Convert each vertex
    }
    return exact_poly;
}

std::vector<Polygon_2> create_and_convert_offset_polygon(double offset_distance, const Polygon_2& polygon_) {
    auto polygon = polygon_;
    if (!polygon.is_counterclockwise_oriented()) {
        polygon.reverse_orientation();
    }

    // Create the offset polygons using Epick kernel
    // create_exterior_skeleton_and_offset_polygons_2()
#if CGAL_VERSION_NR >= 1060000000
#define shared_ptr std::shared_ptr
#else
#define shared_ptr boost::shared_ptr
#endif
    std::vector<shared_ptr<CGAL::Polygon_2<CGAL::Epick>>> offset_polygons;

    if (offset_distance >= 0.) {
        offset_polygons = CGAL::create_exterior_skeleton_and_offset_polygons_2(offset_distance, polygon);
        // erase the first outer frame
        offset_polygons.erase(offset_polygons.begin());
        offset_polygons.front()->reverse_orientation();
    } else {
        offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(-offset_distance, polygon);
    }

    // Convert each offset polygon back to the Epeck kernel
    std::vector<Polygon_2> exact_offset_polygons;
    for (const auto& inexact_poly_ptr : offset_polygons) {
        Polygon_2 exact_poly = convert_polygon(*inexact_poly_ptr);
        exact_offset_polygons.push_back(exact_poly);
    }

    return exact_offset_polygons;
}

template <typename T>
const T& take_first_if_single_item(const std::vector<T>& vec) {
    if (vec.size() == 1) {
        return vec.front();
    }
    throw std::runtime_error("Expected a single item");
}

// Function to write polygons as line segments in OBJ format
void write_polygon_to_obj(std::ofstream& ofs, size_t& vertex_index, bool as_line, const Polygon_2& polygon, const std::string& name) {
    ofs << "o " << name << "\n";  // Object name

    // Write vertices
    for (auto vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit) {
        ofs << "v " << CGAL::to_double(vit->x()) << " " << CGAL::to_double(vit->y()) << " 0\n";
    }

    if (as_line) {
        // Write line segments (edges)
        for (size_t j = 0; j < polygon.size(); ++j) {
            ofs << "l " << vertex_index + j << " " << vertex_index + (j + 1) % polygon.size() << "\n";
        }
    } else {
        ofs << "f";
        for (size_t j = 0; j < polygon.size(); ++j) {
            ofs << " " << vertex_index + j;
        }
        ofs << "\n";
    }

    vertex_index += polygon.size();
}

Polygon_2 circ_to_poly(typename Arrangement_2::Ccb_halfedge_const_circulator circ)
{
    Polygon_2 poly;
    auto curr = circ;
    do {
        poly.push_back(curr->source()->point());
    } while (++curr != circ);
    return poly;
}

Polygon_with_holes_2 circ_to_poly(typename Arrangement_2::Ccb_halfedge_const_circulator circ, typename Arrangement_2::Inner_ccb_const_iterator a, typename Arrangement_2::Inner_ccb_const_iterator b)
{
    Polygon_with_holes_2 poly(circ_to_poly(circ));
    for (auto it = a; it != b; ++it) {
        poly.add_hole(circ_to_poly(*it));
    }
    return poly;
}

void write_polygon_to_svg(std::ostream& ofs, const Polygon_2& polygon) {
    ofs << "<polygon points=\"";
    for (auto vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit) {
        ofs << CGAL::to_double(vit->x()) << "," << CGAL::to_double(vit->y()) << " ";
    }
    ofs << "\" style=\"fill:none;stroke-width:1\" />\n";
}

// Function to write a Polygon_with_holes_2 to an SVG file
void write_polygon_with_holes_to_svg(std::ostream& ofs, const Polygon_with_holes_2& polygon_with_holes) {
    // Write the outer boundary (main polygon)
    if (!polygon_with_holes.is_unbounded()) {
        write_polygon_to_svg(ofs, polygon_with_holes.outer_boundary());
    }

    // Write the holes (if any) with a different color (e.g., red)
    for (auto hit = polygon_with_holes.holes_begin(); hit != polygon_with_holes.holes_end(); ++hit) {
        write_polygon_to_svg(ofs, *hit);
    }
}

void remove_close_points(Polygon_2& p, double eps = 1.e-2) {
    std::vector<CGAL::Point_2<K>> ps;
    ps.reserve(p.size());
    auto I = p.begin();
    auto J = I + 1;
    for (;; ++J) {
        bool last = false;
        if (J == p.end()) {
            J = p.begin();
            last = true;
        }
        // std::cout << "d " << std::sqrt(CGAL::to_double(CGAL::squared_distance(*I, *J))) << std::endl;
        if (CGAL::squared_distance(*I, *J) > (eps * eps)) {
            ps.push_back(*J);
            I = J;
        }
        if (last) {
            break;
        }
    }
    if (ps.size() >= 2 && CGAL::squared_distance(ps.front(), ps.back()) <= eps * eps) {
        // Remove the last point if it is too close to the first point
        ps.pop_back();
    }
    if (ps.size() != p.size()) {
        // std::cerr << "Removed " << (p.size() - ps.size()) << " close points from polygon" << std::endl;
        p = Polygon_2(ps.begin(), ps.end());
    }
}

Polygon_2 fuse_with_offset(const std::vector<Polygon_2>& polygons, double polygon_offset_distance) {
    // Find the outer perimeter using offset - union - negative offset
    std::vector<Polygon_2> offset_polygons;
    for (auto& r : polygons) {
        auto ps = create_and_convert_offset_polygon(polygon_offset_distance, r);
        for (auto& p : ps) {
            if (!p.is_simple()) {
                /*{
                    std::cerr << "[";
                    bool first = true;
                    for (auto& pp : r) {
                        if (!first) {
                            std::cerr << ",";
                        }
                        first = false;
                        std::cerr << "(" << pp.x() << "," << pp.y() << ")";
                    }
                    std::cerr << "]" << std::endl;
                }
                {
                    std::cerr << "[";
                    bool first = true;
                    for (auto& pp : p) {
                        if (!first) {
                            std::cerr << ",";
                        }
                        first = false;
                        std::cerr << "(" << pp.x() << "," << pp.y() << ")";
                    }
                    std::cerr << "]" << std::endl;
                }*/
                throw std::runtime_error("Complex polygon originated from offset");
            }
        }
        offset_polygons.insert(offset_polygons.end(), ps.begin(), ps.end());
    }

    // Perform Boolean union on the offset polygons
    std::vector<Polygon_with_holes_2> unioned_polygons;
    CGAL::join(offset_polygons.begin(), offset_polygons.end(), std::back_inserter(unioned_polygons));
    Polygon_2 fused_removed_close_points;

    {
        std::vector<CGAL::Point_2<K>> ps;
        auto& p = unioned_polygons.front().outer_boundary();
        ps.reserve(p.size());
        auto I = p.begin();
        auto J = I + 1;
        for (;; ++J) {
            bool last = false;
            if (J == p.end()) {
                J = p.begin();
                last = true;
            }
            if (CGAL::squared_distance(*I, *J) > (polygon_offset_distance * polygon_offset_distance)) {
                ps.push_back(*J);
                I = J;
            }
            if (last) {
                break;
            }
        }
        fused_removed_close_points = Polygon_2(ps.begin(), ps.end());
    }

    // Apply negative offset to get the outer perimeter polygon
    auto inner_offset = create_and_convert_offset_polygon(
        // Slightly smaller inset distance for non-manifold situs?
        -polygon_offset_distance + 1.e-8,
        fused_removed_close_points);

    if (inner_offset.size() != 1) {
        throw std::runtime_error("Unexpected union outcome");
    }

    return inner_offset.front();
}

void arrange_cgal_polygons(const std::vector<Polygon_2>& input_polygons_, std::vector<Polygon_2>& output_polygons, double polygon_offset_distance = -1.) {
    if (polygon_offset_distance < 0.) {
        double total_edge_length = 0.;
        size_t num_edges = 0;
        for (auto& p : input_polygons_) {
            for (auto it = p.edges_begin(); it != p.edges_end(); ++it) {
                total_edge_length += std::sqrt(CGAL::to_double(CGAL::squared_distance(it->start(), it->end())));
                num_edges += 1;
            }
        }
        polygon_offset_distance = total_edge_length / num_edges;
    }

    auto input_polygons__ = input_polygons_;
    decltype(input_polygons__) input_polygons;

    for (auto& i : input_polygons__) {
        std::vector<CGAL::Point_2<K>> ps(i.begin(), i.end());
        if (ps.front() == ps.back()) {
            ps.pop_back();
        }
        input_polygons.emplace_back(ps.begin(), ps.end());
    }

    for (auto& polygon : input_polygons) {
        if (!polygon.is_counterclockwise_oriented()) {
            polygon.reverse_orientation();
        }
    }

    for (auto& polygon : input_polygons) {
        remove_close_points(polygon);
    }

#ifdef SVGFILL_DEBUG
    std::ofstream obj("obj.obj");
    size_t vi = 1;

    std::ofstream svg("svg.svg");
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"-10 -10 200 200\">\n";

    for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
        write_polygon_to_obj(obj, vi, true, *it, "input_poly_" + std::to_string(std::distance(input_polygons.begin(), it)));
        write_polygon_to_svg(svg, *it);
    }

    obj << std::flush;
#endif

    typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 2, size_t, CGAL::Box_intersection_d::ID_EXPLICIT> Box;
    std::vector<std::vector<CGAL::Triangle_2<K>>> input_triangulated;

    std::vector<Box> boxes;
    for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
        constexpr double offset = 1.e-3;
        auto b = it->bbox();
        boxes.emplace_back(
            CGAL::Bbox_2(b.xmin() - offset, b.ymin() - offset, b.xmax() + offset, b.ymax() + offset),
            std::distance(input_polygons.begin(), it)
        );

        if (!it->is_simple()) {
#ifdef SVGFILL_DEBUG
            write_polygon_to_obj(obj, vi, true, *it, "self-intersecting");
#endif
            throw std::runtime_error("Self-intersecting input");
        }

        CGAL::Polygon_triangulation_decomposition_2<K> decompositor;
        std::vector<Polygon_2> temp;
        decompositor(*it, std::back_inserter(temp));
        input_triangulated.emplace_back();
        for (auto& pol : temp) {
            auto it = pol.vertices_circulator();
            const auto& p = *(it++);
            const auto& q = *(it++);
            const auto& r = *(it++);
            input_triangulated.back().emplace_back(p, q, r);
        }
    }

    std::set<std::pair<size_t, size_t>> overlaps;

    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), [&input_triangulated, &overlaps](const Box& a, const Box& b) {
        for (auto& t1 : input_triangulated[a.handle()]) {
            bool registered_overlap = false;
            for (auto& t2 : input_triangulated[b.handle()]) {
                if (CGAL::squared_distance(t1, t2) < (1.e-3 * 1.e-3)) {
                    overlaps.insert({
                        (a.handle() < b.handle()) ? a.handle() : b.handle(),
                        (a.handle() < b.handle()) ? b.handle() : a.handle()
                        });
                    registered_overlap = true;
                    break;
                }
            }
            if (registered_overlap) {
                // no need to check other triangles
                break;
            }
        }
    });

    if (true) {
        // solve overlaps by means of subtraction
        // loop over overlaps and subtract the smaller polygon from the larger one

        std::set<size_t> eliminated_polies;
        std::map<size_t, size_t> overlap_counts;
        for (auto& p : overlaps) {
            overlap_counts[p.first]++;
            overlap_counts[p.second]++;
        }

        for (const auto& edge : overlaps) {
            // Skip eliminated
            if (eliminated_polies.find(edge.first) != eliminated_polies.end() ||
                eliminated_polies.find(edge.second) != eliminated_polies.end()) {
                continue;
            }

            // Many overlaps indicate an aggregated polygon, skip them
            /*
            if (overlap_counts[edge.first] > 10 || overlap_counts[edge.second] > 10) {
                if (overlap_counts[edge.first] > 10) {
                    eliminated_polies.insert(edge.first);
                }
                if (overlap_counts[edge.second] > 10) {
                    eliminated_polies.insert(edge.second);
                }
                continue;
            }
            */

            auto& poly1 = input_polygons[edge.first];
            auto& poly2 = input_polygons[edge.second];

            // Skip polygons that have a very high intersection over union
            // ratio, which indicates that they are very likely duplicates
            if (CGAL::do_intersect(poly1, poly2)) {
                std::vector<Polygon_with_holes_2> result;
                CGAL::intersection(poly1, poly2, std::back_inserter(result));
                typename K::FT intersection_area = 0;
                for (auto& r : result) {
                    auto poly_area = r.outer_boundary().area();
                    for (auto& h : r.holes()) {
                        poly_area -= h.area();
                    }
                    intersection_area += poly_area;
                }
                CGAL::Polygon_with_holes_2<K> poly12;
                CGAL::join(poly1, poly2, poly12);
                typename K::FT union_area = poly12.outer_boundary().area();
                for (auto& h : poly12.holes()) {
                    union_area -= h.area();
                }
                if (union_area > 0 && intersection_area / union_area > 0.99) {
                    // std::cerr << intersection_area / union_area << std::endl;
                    eliminated_polies.insert(edge.first);
                    continue;
                }
            }

            if (!(poly1.is_simple() && poly2.is_simple())) {
                continue;
            }

            {
                std::vector<Polygon_with_holes_2> result;
                int assign_to;
                if (poly1.area() > poly2.area()) {
                    CGAL::difference(poly1, take_first_if_single_item(create_and_convert_offset_polygon(1.e-2, poly2)), std::back_inserter(result));
                    poly2 = take_first_if_single_item(create_and_convert_offset_polygon(-1.e-2, poly2));
                    assign_to = edge.first;
                } else {
                    CGAL::difference(poly2, take_first_if_single_item(create_and_convert_offset_polygon(1.e-2, poly1)), std::back_inserter(result));
                    poly1 = take_first_if_single_item(create_and_convert_offset_polygon(-1.e-2, poly1));
                    assign_to = edge.second;
                }
                if (result.size() == 1) {
                    input_polygons[assign_to] = result.front().outer_boundary();
                    if (result.front().number_of_holes() == 0) {

                    } else {
                        /*
                        write_polygon_to_obj(obj, vi, true, result.front().outer_boundary(), "invalid_outer");
                        size_t ii = 0;
                        for (auto& i : result.front().holes()) {
                            write_polygon_to_obj(obj, vi, true, i, "invalid_hole_" + std::to_string(ii + 1));
                        }
                        obj << std::flush;
                        */
                        eliminated_polies.insert(assign_to == edge.first ? edge.second : edge.first);
                    }
                } else {
#ifdef SVGFILL_DEBUG
                    write_polygon_to_obj(obj, vi, true, poly1, "err1");
                    write_polygon_to_obj(obj, vi, true, poly2, "err2");
                    obj << std::flush;
#endif
                    throw std::runtime_error("Unexpected union outcome");
                }
            }
        }

        // iterate over the eliminated polygons and remove them from the input polygons
        for (auto it = eliminated_polies.rbegin(); it != eliminated_polies.rend(); ++it) {
            input_polygons.erase(input_polygons.begin() + *it);
        }
    }

    if (false) {
        // solve overlap by means of union into components
        std::vector<std::vector<size_t>> adj(input_polygons.size());
        for (const auto& edge : overlaps) {
            adj[edge.first].push_back(edge.second);
            adj[edge.second].push_back(edge.first);
        }

        std::vector<bool> visited(input_polygons.size(), false);
        std::vector<std::vector<size_t>> connected_components;

        for (size_t v = 0; v < input_polygons.size(); ++v) {
            if (!visited[v]) {
                connected_components.emplace_back();

                std::stack<size_t> stack;
                stack.push(v);
                visited[v] = true;

                while (!stack.empty()) {
                    size_t u = stack.top();
                    stack.pop();
                    connected_components.back().push_back(u);

                    for (size_t neighbor : adj[u]) {
                        if (!visited[neighbor]) {
                            visited[neighbor] = true;
                            stack.push(neighbor);
                        }
                    }
                }
            }
        }

        std::vector<Polygon_2> fused_polies;

        for (auto& comp : connected_components) {
            std::vector<Polygon_2> comp_polies;
            if (comp.size() == 1) {
                fused_polies.push_back(input_polygons[comp.front()]);
            } else {
                for (auto& c : comp) {
                    comp_polies.push_back(input_polygons[c]);
                }
                fused_polies.push_back(fuse_with_offset(comp_polies, 1.e-2));
            }
        }

#ifdef SVGFILL_DEBUG
        for (auto it = fused_polies.begin(); it != fused_polies.end(); ++it) {
            write_polygon_to_obj(obj, vi, true, *it, "fused_poly_" + std::to_string(std::distance(fused_polies.begin(), it)));
            write_polygon_to_svg(svg, *it);
        }
#endif

        input_polygons = fused_polies;
    }

#ifdef SVGFILL_DEBUG
    for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
        write_polygon_to_obj(obj, vi, true, *it, "processed_input_poly_" + std::to_string(std::distance(input_polygons.begin(), it)));
        write_polygon_to_svg(svg, *it);
    }
#endif

    auto input_polygon_boundary = [&input_polygons](const CGAL::Point_2<K>& p, double tol = 1.e-5) {
        // unfortunately some imprecision slept into the code so we can't
        // so we can't just use has_on_boundary() anymore
        double D = std::numeric_limits<double>::infinity();
        for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
            for (auto jt = it->edges_begin(); jt != it->edges_end(); ++jt) {
                const auto& seg = *jt;
                auto d = std::sqrt(CGAL::to_double(CGAL::squared_distance(seg, p)));
                if (d < D) {
                    D = d;
                }
                if (d < tol) {
                    return it;
                }
            }
        }
        return input_polygons.end();
    };

    auto close_input_point = [&input_polygons](const CGAL::Point_2<K>& P) {
        CGAL::Point_2<K> closest;
        double closest_distance = std::numeric_limits<double>::infinity();
        auto input_it = input_polygons.end();

        // unfortunately some imprecision slept into the code so we can't
        // so we can't just use has_on_boundary() anymore
        for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
            for (auto& p : *it) {
                auto d = std::sqrt(CGAL::to_double(CGAL::squared_distance(P, p)));
                if (d < closest_distance) {
                    closest_distance = d;
                    closest = p;
                    input_it = it;
                }
            }
        }

        return std::make_pair(input_it, closest);
    };

    auto project_input_point = [&input_polygons](const CGAL::Point_2<K>& P) {
        CGAL::Point_2<K> closest;
        typename K::FT closest_sq_distance = std::numeric_limits<double>::infinity();
        auto input_it = input_polygons.end();

        // unfortunately some imprecision slept into the code so we can't
        // so we can't just use has_on_boundary() anymore
        for (auto it = input_polygons.begin(); it != input_polygons.end(); ++it) {
            for (auto jt = it->edges_begin(); jt != it->edges_end(); ++jt) {
                auto Pp = jt->supporting_line().projection(P);
                auto d = CGAL::squared_distance(Pp, P);
                if (d < closest_sq_distance) {
                    closest_sq_distance = d;
                    closest = Pp;
                    input_it = it;
                }
            }
        }

        return std::make_pair(input_it, closest);
    };

    // Find the outer perimeter using offset - union - negative offset
    std::vector<Polygon_2> offset_polygons;
    for (auto& r : input_polygons) {
        auto R = r;
        if (!R.is_counterclockwise_oriented()) {
            R.reverse_orientation();
        }

        // Overlap removal can also result in close points causing problems when converted into non-exact nt
        remove_close_points(R);

        auto ps = create_and_convert_offset_polygon(polygon_offset_distance, R);
        for (auto& p : ps) {
            if (!p.is_simple()) {
                /*{
                    std::cerr << "input [";
                    bool first = true;
                    for (auto& pp : r) {
                        if (!first) {
                            std::cerr << ",";
                        }
                        first = false;
                        std::cerr << "(" << pp.x() << "," << pp.y() << ")";
                    }
                    std::cerr << "]" << std::endl;
                }

                {
                    std::cerr << "[";
                    bool first = true;
                    for (auto& pp : p) {
                        if (!first) {
                            std::cerr << ",";
                        }
                        first = false;
                        std::cerr << "(" << pp.x() << "," << pp.y() << ")";
                    }
                    std::cerr << "]" << std::endl;
                }*/

                throw std::runtime_error("Complex polygon originated from offset");
            }
        }
        offset_polygons.insert(offset_polygons.end(), ps.begin(), ps.end());
    }

#ifdef SVGFILL_DEBUG
    for (auto it = offset_polygons.begin(); it != offset_polygons.end(); ++it) {
        write_polygon_to_obj(obj, vi, true, *it, "offset_poly_" + std::to_string(std::distance(offset_polygons.begin(), it)));
        write_polygon_to_svg(svg, *it);
    }
#endif

    // Perform Boolean union on the offset polygons
    std::vector<Polygon_with_holes_2> unioned_polygons;
    CGAL::join(offset_polygons.begin(), offset_polygons.end(), std::back_inserter(unioned_polygons));

#ifdef SVGFILL_DEBUG
    write_polygon_to_obj(obj, vi, true, unioned_polygons.front().outer_boundary(), "offset_poly_joined");
    write_polygon_to_svg(svg, unioned_polygons.front().outer_boundary());

#endif

    Polygon_2 fused_removed_close_points;
    {
        std::vector<CGAL::Point_2<K>> ps;
        auto& p = unioned_polygons.front().outer_boundary();
        ps.reserve(p.size());
        auto I = p.begin();
        auto J = I + 1;
        for (;; ++J) {
            bool last = false;
            if (J == p.end()) {
                J = p.begin();
                last = true;
            }
            // if (CGAL::squared_distance(*I, *J) > (polygon_offset_distance * polygon_offset_distance)) {
            if (CGAL::squared_distance(*I, *J) > (1.e-4 * 1.e-4)) {
                ps.push_back(*J);
                I = J;
            }
            if (last) {
                break;
            }
        }
        fused_removed_close_points = Polygon_2(ps.begin(), ps.end());
    }

    // Apply negative offset to get the outer perimeter polygon
    auto inner_offset = create_and_convert_offset_polygon(
        // Because polygon_offset is inexact, make sure our inset distance is slightly larger
        // std::nexttoward(-polygon_offset_distance, -std::numeric_limits<double>::infinity()),

        // 1.e-8 even was too large and still resulted in slivers of triangle around the perimeter  
        -polygon_offset_distance - 1.e-5,
        fused_removed_close_points);

#ifdef SVGFILL_DEBUG
    write_polygon_to_obj(obj, vi, true, inner_offset.front(), "joined_inset");
    write_polygon_to_svg(svg, inner_offset.front());
#endif

    // there is non-insignificant chance that around the outer boundary, vertices are located in
    // between of the input polyhedra, but intermediate vertices result in triangles that will no longer
    // span between the two spaces with two edges and therefore cause the topological centre line
    // to no run up to the center. Eliminate all vertices that are not on the polyhedral boundary of polygon.

    // this theory proved to be false. once we have topological end points in our graph that are
    // connected to input polyhedra to form closed cells, we move those topological end points to
    // the average of the input polyhedra corner points, thus effectively also moving them outwards.
    {
        for (auto& i : inner_offset) {
            std::vector<CGAL::Point_2<K>> ps;
            for (auto& p : i) {
                if (input_polygon_boundary(p, 1.e-3) != input_polygons.end()) {
                    ps.push_back(p);
                }
            }
            i = Polygon_2(ps.begin(), ps.end());
        }
    }

#ifdef SVGFILL_DEBUG
    write_polygon_to_obj(obj, vi, true, inner_offset.front(), "joined_inset_cleaned");
    write_polygon_to_svg(svg, inner_offset.front());
#endif

    // Subtract original polygons from outer perimeter
    std::vector<Polygon_with_holes_2> difference_result, difference_result_subdivided;
    for (auto& i : inner_offset) {
        std::vector<Polygon_with_holes_2> working_copy;
        working_copy.emplace_back(i);

        for (auto& r : input_polygons) {
            std::vector<Polygon_with_holes_2> temp_working_copy;
            for (auto& wc : working_copy) {
                CGAL::difference(wc, r, std::back_inserter(temp_working_copy));
            }
            working_copy = temp_working_copy;
        }
        difference_result.insert(difference_result.end(), working_copy.begin(), working_copy.end());
    }

    // subdivide difference_result to have better behave triangulation

    {
        const double max_distance = polygon_offset_distance / 8.;
        auto subdivide_polygon = [max_distance](const Polygon_2& p) {
            std::vector<Point_2> points;
            for (auto it = p.edges_begin(); it != p.edges_end(); ++it) {
                const auto& seg = *it;
                auto num_splits = (int)std::ceil(std::sqrt(CGAL::to_double(seg.squared_length())) / max_distance) - 1;
                points.push_back(seg.source());
                for (auto i = 0; i < num_splits; ++i) {
                    auto d = (seg.target() - seg.source()) / (num_splits + 1) * (i + 1);
                    points.push_back(seg.source() + d);
                }
            }
            return Polygon_2(points.begin(), points.end());
        };

        for (auto& pwh : difference_result) {
            // Subdivide outer boundary
            Polygon_2 outer = subdivide_polygon(pwh.outer_boundary());
            // Subdivide holes
            std::vector<Polygon_2> holes;
            for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
                holes.push_back(subdivide_polygon(*hit));
            }
            // Construct new Polygon_with_holes_2
            difference_result_subdivided.push_back(Polygon_with_holes_2(outer, holes.begin(), holes.end()));
        }
    }

#ifdef SVGFILL_DEBUG
    for (auto it = difference_result_subdivided.begin(); it != difference_result_subdivided.end(); ++it) {
        auto i = std::distance(difference_result_subdivided.begin(), it);
        write_polygon_to_obj(obj, vi, true, it->outer_boundary(), "difference_result_subdivided_" + std::to_string(i));
        write_polygon_to_svg(svg, it->outer_boundary());
        for (auto& p : it->holes()) {
            write_polygon_to_obj(obj, vi, true, p, "difference_result_subdivided_" + std::to_string(i));
            write_polygon_to_svg(svg, p);
        }
    }
#endif

    std::list<CGAL::Polygon_2<K>> triangular_polygons;

    for (auto& pwh : difference_result_subdivided) {
        CGAL::Polygon_triangulation_decomposition_2<K> decompositor;
        decompositor(pwh, std::back_inserter(triangular_polygons));
    }

    triangular_polygons.erase(std::remove_if(triangular_polygons.begin(), triangular_polygons.end(), [](const CGAL::Polygon_2<K>& p) {
        return CGAL::to_double(p.area()) < 1.e-8;
    }), triangular_polygons.end());

#ifdef SVGFILL_DEBUG
    for (auto it = triangular_polygons.begin(); it != triangular_polygons.end(); ++it) {
        write_polygon_to_obj(obj, vi, false, *it, "tri_" + std::to_string(std::distance(triangular_polygons.begin(), it)));
        write_polygon_to_svg(svg, *it);
    }
#endif

    // Build maps of triangle -> edge and edge -> triangle in order to do traversal on the 'corridor mesh'
    std::map<std::pair<Point_2, Point_2>, std::vector<CGAL::Polygon_2<K>*>> segment_to_facet;
    std::map<std::pair<Point_2, Point_2>, Point_2> segment_to_midpoint;
    std::map<Point_2, std::pair<Point_2, Point_2>> midpoint_to_segment;
    std::map<CGAL::Polygon_2<K>*, std::vector<std::pair<Point_2, Point_2>>> facet_to_segment;

    for (auto& tri : triangular_polygons) {
        for (size_t i = 0; i < 3; ++i) {
            size_t j = (i + 1) % 3;
            auto& pi = tri.vertex(i);
            auto& pj = tri.vertex(j);
            const bool orientation = std::lexicographical_compare(pi.cartesian_begin(), pi.cartesian_end(), pj.cartesian_begin(), pj.cartesian_end());
            std::pair<Point_2, Point_2> seg(orientation ? pi : pj, orientation ? pj : pi);
            segment_to_facet[seg].push_back(&tri);
            facet_to_segment[&tri].push_back(seg);
        }
    }

    // Register midpoints on the edges within the 'corridor mesh' that span multiple input polygons
    for (auto& p : segment_to_facet) {
        auto center = CGAL::ORIGIN + (((p.first.first - CGAL::ORIGIN) + (p.first.second - CGAL::ORIGIN)) / 2);


        auto p1index = input_polygon_boundary(p.first.first);
        auto p2index = input_polygon_boundary(p.first.second);
        if (p1index != input_polygons.end() && p2index != input_polygons.end() && p1index != p2index) {
            segment_to_midpoint[p.first] = center;
            midpoint_to_segment[center] = p.first;
        }

        if (p1index != input_polygons.end() && p2index != input_polygons.end() && p1index != p2index) {
            segment_to_midpoint[p.first] = center;
            midpoint_to_segment[center] = p.first;
        }
    }

#ifdef SVGFILL_DEBUG
    obj << "o network_1\n";
#endif

    // Observe corridor mesh topology to join edge midpoints into a network
    std::map<Point_2, std::vector<Point_2>> line_graph;
    for (auto& p : segment_to_midpoint) {
        for (auto& q : segment_to_facet[p.first]) {
            for (auto& r : facet_to_segment[q]) {
                if (p.first == r) {
                    continue;
                }
                decltype(segment_to_midpoint)::const_iterator it;
                if ((it = segment_to_midpoint.find(r)) != segment_to_midpoint.end()) {
                    line_graph[p.second].push_back(it->second);

#ifdef SVGFILL_DEBUG
                    obj << "v " << CGAL::to_double(p.second.x()) << " " << CGAL::to_double(p.second.y()) << " 0\n";
                    obj << "v " << CGAL::to_double(it->second.x()) << " " << CGAL::to_double(it->second.y()) << " 0\n";
                    obj << "l " << vi++;
                    obj << " " << vi++ << "\n";

                    svg << "<line x1=\"" << CGAL::to_double(p.second.x()) << "\" y1=\"" << CGAL::to_double(p.second.y()) << "\" x2=\"" << CGAL::to_double(it->second.x()) << "\" y2=\"" << CGAL::to_double(it->second.y()) << "\" />";
#endif
                }
            }
        }
    }

    // Find triangles in this network often occuring at junctions in the corridor mesh
    std::set<Triangle<K>> triangles;
    std::function<void(std::vector<Point_2>&)> find_triangles_recursive;
    find_triangles_recursive = [&](std::vector<Point_2>& path) -> void {
        // If depth reaches 3, check for a triangle
        if (path.size() == 3) {
            // Check if we can complete the triangle by going from the current point back to the start
            const std::vector<Point_2>& neighbors_current = line_graph.at(path.back());
            if (std::find(neighbors_current.begin(), neighbors_current.end(), path.front()) != neighbors_current.end()) {
                // We found a triangle, add it to the set
                Triangle<K> triangle = { path[0], path[1], path[2] };
                std::sort(triangle.begin(), triangle.end());
                triangles.insert(triangle);
            }
            return;
        }

        // Otherwise, continue exploring neighbors
        const std::vector<Point_2>& neighbors = line_graph.at(path.back());
        for (const Point_2& neighbor : neighbors) {
            if (std::find(path.begin(), path.end(), neighbor) == path.end()) {
                path.push_back(neighbor);
                find_triangles_recursive(path);
                path.pop_back();  // Backtrack
            }
        }
    };

    for (auto& p : line_graph) {
        std::vector<Point_2> ps = { p.first };
        find_triangles_recursive(ps);
    }

    // For every triangle found in the network we eliminate one edge to break the cycle
    // The edge we eliminate is the edge with the greatest angle with any of it's neighbours

    // non exact time, we need sqrt
    using SK = CGAL::Simple_cartesian<double>;
    CGAL::Cartesian_converter<K, SK> C{};

#ifdef SVGFILL_DEBUG
    obj << "o eliminated\n";
#endif

    std::set<std::pair<Point_2, Point_2>> eliminated_segments;

    for (auto& t : triangles) {
        Triangle<SK> st;
        std::transform(t.begin(), t.end(), st.begin(), C);

        double global_min_abs_dot = std::numeric_limits<double>::infinity();
        size_t global_min_abs_dot_index;

        for (size_t i = 0; i < 3; ++i) {
            auto j = (i + 2) % 3;
            auto e0 = st[i] - st[j];
            e0 /= std::sqrt(e0.squared_length());

            double max_abs_dot = 0.;

            {
                auto& ni = line_graph[t[i]];
                for (auto& n : ni) {
                    if (std::find(t.begin(), t.end(), n) == t.end()) {
                        // not contained in triangle
                        auto sn = C(n);
                        auto en = sn - st[i];
                        en /= std::sqrt(en.squared_length());
                        auto dot = std::abs(en * e0);

                        if (dot > max_abs_dot) {
                            max_abs_dot = dot;
                        }
                    }
                }
            }

            {
                auto& nj = line_graph[t[j]];
                for (auto& n : nj) {
                    if (std::find(t.begin(), t.end(), n) == t.end()) {
                        // not contained in triangle
                        auto sn = C(n);
                        auto en = sn - st[j];
                        en /= std::sqrt(en.squared_length());
                        auto dot = std::abs(en * e0);

                        if (dot > max_abs_dot) {
                            max_abs_dot = dot;
                        }
                    }
                }

            }

            if (max_abs_dot < global_min_abs_dot) {
                global_min_abs_dot = max_abs_dot;
                global_min_abs_dot_index = i;
            }
        }

        {
            auto i = global_min_abs_dot_index;
            auto j = (i + 2) % 3;

            eliminated_segments.insert({ t[i], t[j] });
            eliminated_segments.insert({ t[j], t[i] });

#ifdef SVGFILL_DEBUG
            obj << "v " << st[j].x() << " " << st[j].y() << " 0\n";
            obj << "v " << st[i].x() << " " << st[i].y() << " 0\n";
            obj << "l " << vi++;
            obj << " " << vi++ << "\n";

            svg << "<line x1=\"" << CGAL::to_double(st[j].x()) << "\" y1=\"" << CGAL::to_double(st[j].y()) << "\" x2=\"" << CGAL::to_double(st[i].x()) << "\" y2=\"" << CGAL::to_double(st[i].y()) << "\" stroke=\"red\" />";
#endif
        }

    }

    Graph2D<K> G2(line_graph);
    for (auto& e : eliminated_segments) {
        G2.remove_edge(e.first, e.second);
    }

    auto G = G2.weld_vertices();

#ifdef SVGFILL_DEBUG
    obj << "o network_2\n";
    for (auto it = G.edges_begin(); it != G.edges_end(); ++it) {
        obj << "v " << CGAL::to_double(it->first.x()) << " " << CGAL::to_double(it->first.y()) << " 0\n";
        obj << "v " << CGAL::to_double(it->second.x()) << " " << CGAL::to_double(it->second.y()) << " 0\n";
        obj << "l " << vi++;
        obj << " " << vi++ << "\n";
    }
    obj << std::flush;
#endif

    auto is_parallel_2degree_node = [](decltype(G)::vertex_const_iterator vit) {
        auto it = vit->second.begin();
        auto& P = *it++;
        auto& Q = *it++;
        auto e1 = P - vit->first;
        auto e2 = vit->first - Q;
        if (e1.squared_length() == 0 || e2.squared_length() == 0) {
            // @todo why does this happen?
            return false;
        }
        e1 /= std::sqrt(CGAL::to_double(e1.squared_length()));
        e2 /= std::sqrt(CGAL::to_double(e2.squared_length()));
        return std::abs(CGAL::to_double(e1 * e2)) > (1. - 1.e-5);
    };

    {
        // Remove colinear vertices
        size_t n_vertices_removed = 0;
        for (auto vit = G.vertices_begin(); vit != G.vertices_end();) {
            if (vit->second.size() == 2) {
                if (is_parallel_2degree_node(vit)) {
                    vit = G.eliminate_vertex(vit);
                    ++n_vertices_removed;
                } else {
                    ++vit;
                }
            } else {
                ++vit;
            }
        }
        // std::cout << "Eliminated " << n_vertices_removed << " vertices" << std::endl;
    }

    // Ortho edge slide
    {
        std::list<CGAL::Segment_2<K>> edges_to_remove, edges_to_insert;

        for (auto vit = G.vertices_begin(); vit != G.vertices_end(); ++vit) {
            auto& selected = vit->first;

            if (vit->second.size() >= 3) {
                for (auto vjt = vit->second.begin(); vjt != vit->second.end(); ++vjt) {
                    auto& neighbour = *vjt;
                    bool processed_neighbour = false;

                    if (G.find(neighbour)->second.size() == 2 && !is_parallel_2degree_node(G.find(neighbour))) {
                        auto vkt = G.find(neighbour)->second.begin();
                        if (selected == *vkt) {
                            vkt++;
                        }
                        auto& other = *vkt;

                        if ((other - neighbour).squared_length() < (neighbour - selected).squared_length()) {
                            continue;
                        }

                        auto incoming = CGAL::Ray_2<K>(other, neighbour - other);
                        boost::optional<CGAL::Segment_2<K>> closest_neighbouring_segment;
                        boost::optional<CGAL::Point_2<K>> closest_intersection_point;
                        K::FT sq_distance_along_ray = std::numeric_limits<double>::infinity();

                        for (auto vlt = vit->second.begin(); vlt != vit->second.end(); ++vlt) {
                            auto& other_neighbour = *vlt;
                            if (vlt != vjt) {
                                CGAL::Segment_2<K> neighbouring_segment(selected, other_neighbour);
                                auto x = CGAL::intersection(incoming, neighbouring_segment);
                                if (x) {
#if CGAL_VERSION_NR >= 1060000000
#define variant_get std::get_if
#else
#define variant_get boost::get
#endif

                                    if (auto* xp = variant_get<CGAL::Point_2<K>>(&*x)) {
                                        auto dist = ((*xp) - other).squared_length();
                                        if (dist < sq_distance_along_ray) {
                                            closest_neighbouring_segment = neighbouring_segment;
                                            closest_intersection_point = *xp;
                                            sq_distance_along_ray = dist;
                                        }
                                    }
                                }
                            }
                        }

                        if (closest_intersection_point && closest_neighbouring_segment) {
                            edges_to_remove.push_back(*closest_neighbouring_segment);
                            edges_to_remove.push_back({ neighbour, selected });
                            edges_to_insert.push_back({ closest_neighbouring_segment->source(), *closest_intersection_point });
                            edges_to_insert.push_back({ closest_neighbouring_segment->target(), *closest_intersection_point });
                            edges_to_insert.push_back({ neighbour, *closest_intersection_point });

                            processed_neighbour = true;
                        }
                    }
                    if (processed_neighbour) {
                        // Only one neigbour is processed because otherwise we obtain intersections
                        break;
                    }
                }
            }
        }

        for (auto& s : edges_to_remove) {
            G.remove_edge(s.source(), s.target());
        }


        for (auto& s : edges_to_insert) {
            G.insert(s.source(), s.target());
        }

#ifdef SVGFILL_DEBUG
        obj << "o network_3\n";
        for (auto it = G.edges_begin(); it != G.edges_end(); ++it) {
            obj << "v " << CGAL::to_double(it->first.x()) << " " << CGAL::to_double(it->first.y()) << " 0\n";
            obj << "v " << CGAL::to_double(it->second.x()) << " " << CGAL::to_double(it->second.y()) << " 0\n";
            obj << "l " << vi++;
            obj << " " << vi++ << "\n";
        }
#endif
    }

    // Now plot the edges on an arrangement in order to find planar cycles
    // and merge the corridor-halves with their neighbouring input polygon

    Arrangement_2 arr;

    for (auto it = G.edges_begin(); it != G.edges_end(); ++it) {
        CGAL::insert(arr, Segment_2(it->first, it->second));
    }

    std::list<std::pair<Point_2, Point_2>> move_ops;
    std::list<std::pair<Point_2, Point_2>> edge_ops;

    for (auto it = G.vertices_begin(); it != G.vertices_end(); ++it) {
        if (it->second.size() == 1) {
            auto& M = it->first;

            decltype(midpoint_to_segment)::mapped_type* q = nullptr;

            if (midpoint_to_segment.find(M) == midpoint_to_segment.end()) {
                typename K::FT min_sq_distance = std::numeric_limits<double>::infinity();
                for (auto& pa : midpoint_to_segment) {
                    if (CGAL::squared_distance(pa.first, M) < min_sq_distance) {
                        q = &pa.second;
                        min_sq_distance = CGAL::squared_distance(pa.first, M);
                    }
                }
            } else {
                q = &midpoint_to_segment[M];
            }

            if (q == nullptr) {
                continue;
            }

            // distance from unioned - shoot ray?

            // for now we choose to map point to the midpoint of the found two close points.

            auto pq = close_input_point(q->first);
            auto pr = close_input_point(q->second);

            auto Q = pq.second;
            auto R = pr.second;

            if (Q == R) {
                // this can happen in situations like this:
                // where Q and R are co-located, because the point R' is further away
                // in that case M + M-Q should gives is x that we then project onto the
                // input boundary
                // 
                // 
                // ┌───────┐                
                // │       │                
                // │       │                
                // │       │                
                // └───────o   <--Q,R                
                //                          
                // ────────o   <--M             
                //                          
                // ┌───────x───────────────o  <---R'
                // │                       │
                // │                       │
                // │                       │
                // │                       │
                // └───────────────────────┘

                // @todo is this projection actually necessary or is it already 'exact enough'?
                R = project_input_point(M + (M - Q)).second;
            }

            auto avg = CGAL::ORIGIN + ((Q - CGAL::ORIGIN) + (R - CGAL::ORIGIN)) / 2;

            move_ops.push_front({ M, avg });
            edge_ops.push_front({ avg, Q });
            edge_ops.push_front({ avg, R });
        }
    }

#ifdef SVGFILL_DEBUG
    obj << "o network_4\n";
#endif

    // note that we actually don't move but draw an edge
    for (auto& pq : move_ops) {
        CGAL::insert(arr, Segment_2(pq.first, pq.second));

#ifdef SVGFILL_DEBUG
        obj << "v " << CGAL::to_double(pq.first.x()) << " " << CGAL::to_double(pq.first.y()) << " 0\n";
        obj << "v " << CGAL::to_double(pq.second.x()) << " " << CGAL::to_double(pq.second.y()) << " 0\n";
        obj << "l " << vi++;
        obj << " " << vi++ << "\n";
#endif
    }


    for (auto& pq : edge_ops) {
        CGAL::insert(arr, Segment_2(pq.first, pq.second));

#ifdef SVGFILL_DEBUG
        obj << "v " << CGAL::to_double(pq.first.x()) << " " << CGAL::to_double(pq.first.y()) << " 0\n";
        obj << "v " << CGAL::to_double(pq.second.x()) << " " << CGAL::to_double(pq.second.y()) << " 0\n";
        obj << "l " << vi++;
        obj << " " << vi++ << "\n";
#endif
    }

    // Plot input polygons
    for (auto& poly : input_polygons) {
        for (size_t i = 0; i != poly.size(); ++i) {
            auto j = (i + 1) % poly.size();
            CGAL::insert(arr, Segment_2(poly.vertex(i), poly.vertex(j)));
        }
    }


#ifdef SVGFILL_DEBUG
    {
        obj << "o arrangement_1\n";
        for (auto it = arr.edges_begin(); it != arr.edges_end(); ++it) {
            auto& p = it->source()->point();
            auto& q = it->target()->point();
            obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
            obj << "v " << CGAL::to_double(q.x()) << " " << CGAL::to_double(q.y()) << " 0\n";
            obj << "l " << vi++;
            obj << " " << vi++ << "\n";
        }
    }
#endif

    /* {
        // debug, add outer bounds so that we can plot the face for any remaining edges
        auto poly = unioned_polygons.front().outer_boundary();
        for (size_t i = 0; i != poly.size(); ++i) {
            auto j = (i + 1) % poly.size();
            CGAL::insert(arr, Segment_2(poly.vertex(i), poly.vertex(j)));
        }
    } */

    // Now loop over the arrangement faces, when a face coincides with a point on the
    // corridor network we know it needs to be joined with an input polygon. In that
    // case the edges need to be eliminated that correspond to original geometry.

    size_t face_id = 0;
    std::set<Arrangement_2::Halfedge_handle> edges_to_remove;

    for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) {
        if (it->is_unbounded()) {
            continue;
        }
        bool is_corridor = false;
        {
            auto curr = it->outer_ccb();
            do {
                auto& p = curr++->source()->point();
                if (G.find(p) != G.vertices_end()) {
                    is_corridor = true;
                    break;
                }
            } while (curr != it->outer_ccb());

            for (auto jt = it->inner_ccbs_begin(); jt != it->inner_ccbs_end(); ++jt) {
                curr = *jt;
                do {
                    auto& p = curr++->source()->point();
                    if (G.find(p) != G.vertices_end()) {
                        is_corridor = true;
                        break;
                    }
                } while (curr != *jt);
                if (is_corridor) {
                    break;
                }
            }
        }

        if (is_corridor) {
            auto curr = it->outer_ccb();
            do {
                auto& p = curr->source()->point();
                auto& q = curr->target()->point();
                auto center = CGAL::ORIGIN + (((p - CGAL::ORIGIN) + (q - CGAL::ORIGIN)) / 2);
                auto p1index = input_polygon_boundary(center);
                const bool on_orig_bound = p1index != input_polygons.end();
                if (on_orig_bound) {
                    if (edges_to_remove.find(curr->twin()) != edges_to_remove.end()) {
                        // std::cerr << "Warning trying to delete edge twice" << std::endl;
                    } else {
                        edges_to_remove.insert(curr);
                    }
                }
                curr++;
            } while (curr != it->outer_ccb());

            for (auto jt = it->inner_ccbs_begin(); jt != it->inner_ccbs_end(); ++jt) {
                curr = *jt;
                do {
                    auto& p = curr->source()->point();
                    auto& q = curr->target()->point();
                    auto center = CGAL::ORIGIN + (((p - CGAL::ORIGIN) + (q - CGAL::ORIGIN)) / 2);
                    auto p1index = input_polygon_boundary(center);
                    const bool on_orig_bound = p1index != input_polygons.end();
                    if (on_orig_bound) {
                        if (edges_to_remove.find(curr->twin()) != edges_to_remove.end()) {
                            // std::cerr << "Warning trying to delete edge twice" << std::endl;
                        } else {
                            edges_to_remove.insert(curr);
                        }
                    }
                    curr++;
                } while (curr != *jt);
                if (is_corridor) {
                    break;
                }
            }
        }

#ifdef SVGFILL_DEBUG
        write_polygon_to_svg(svg, circ_to_poly(it->outer_ccb()));

        obj << "o " << "face_";
        if (is_corridor) {
            obj << "corri_";
        }
        obj << face_id++ << "\n";

        std::ostringstream oss;

        {
            auto vv = vi;
            auto curr = it->outer_ccb();
            do {
                auto& p = curr->source()->point();
                obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
                oss << "l " << vi++;
                ++curr;
                if (curr == it->outer_ccb()) {
                    oss << " " << vv << "\n";
                } else {
                    oss << " " << vi << "\n";
                }
            } while (curr != it->outer_ccb());
        }

        for (auto jt = it->inner_ccbs_begin(); jt != it->inner_ccbs_end(); ++jt) {
            auto vv = vi;
            auto curr = *jt;
            do {
                auto& p = curr->source()->point();
                obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
                oss << "l " << vi++;
                ++curr;
                if (curr == *jt) {
                    oss << " " << vv << "\n";
                } else {
                    oss << " " << vi << "\n";
                }
            } while (curr != *jt);
        }

        obj << oss.str();
#endif
    }

    size_t remove_id = 0;
    for (auto& e : edges_to_remove) {
#ifdef SVGFILL_DEBUG
        obj << "o " << "remove_" << remove_id++ << "\n";
        {
            auto& p = e->source()->point();
            obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
        }
        {
            auto& p = e->target()->point();
            obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
        }
        obj << "l " << vi++;
        obj << " " << vi++ << std::endl;
#endif
        CGAL::remove_edge(arr, e);
    }

    for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) {
        if (it->is_unbounded()) {
            continue;
        }

        output_polygons.push_back(circ_to_poly(it->outer_ccb()));

#ifdef SVGFILL_DEBUG
        write_polygon_to_svg(svg, circ_to_poly(it->outer_ccb()));

        obj << "o " << "merged_face_";
        obj << face_id++ << "\n";

        std::ostringstream oss;

        {
            auto vv = vi;
            auto curr = it->outer_ccb();
            do {
                auto& p = curr->source()->point();
                obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
                oss << "l " << vi++;
                ++curr;
                if (curr == it->outer_ccb()) {
                    oss << " " << vv << "\n";
                } else {
                    oss << " " << vi << "\n";
                }
            } while (curr != it->outer_ccb());
        }

        for (auto jt = it->inner_ccbs_begin(); jt != it->inner_ccbs_end(); ++jt) {
            auto vv = vi;
            auto curr = *jt;
            do {
                auto& p = curr->source()->point();
                obj << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " 0\n";
                oss << "l " << vi++;
                ++curr;
                if (curr == *jt) {
                    oss << " " << vv << "\n";
                } else {
                    oss << " " << vi << "\n";
                }
            } while (curr != *jt);
        }

        obj << oss.str();
#endif
    }

#ifdef SVGFILL_DEBUG
    svg << "</svg>\n";
#endif
}

#ifndef SVGFILL_MAIN

bool svgfill::arrange_polygons(const std::vector<svgfill::polygon_2>& polygons, std::vector<svgfill::polygon_2>& arranged)
{
    std::vector<Polygon_2> cgal_polygons, cgal_polygons_out;
    std::transform(polygons.begin(), polygons.end(), std::back_inserter(cgal_polygons), [](auto& poly) {
        Polygon_2 result;
        std::transform(poly.boundary.begin(), poly.boundary.end(), std::back_inserter(result), [](auto& p) {
            return Point_2(p[0], p[1]);
        });
        return result;
    });
    arrange_cgal_polygons(cgal_polygons, cgal_polygons_out);
    std::transform(cgal_polygons_out.begin(), cgal_polygons_out.end(), std::back_inserter(arranged), [](auto& poly) {
        svgfill::polygon_2 result;
        std::transform(poly.begin(), poly.end(), std::back_inserter(result.boundary), [](auto& pt) {
            return svgfill::point_2{
                CGAL::to_double(pt.cartesian(0)),
                CGAL::to_double(pt.cartesian(1)),
            };
        });
        return result;
    });
    return true;
}

#else

template <typename T>
Polygon_2 create_rectangle(T x_min, T y_min, T x_max, T y_max) {
    Polygon_2 rectangle;
    rectangle.push_back(Point_2(x_min, y_min));
    rectangle.push_back(Point_2(x_max, y_min));
    rectangle.push_back(Point_2(x_max, y_max));
    rectangle.push_back(Point_2(x_min, y_max));
    return rectangle;
}

#include <nlohmann/json.hpp>

int main(int argc, char** argv) {
    std::vector<Polygon_2> input_polygons, output;

    if (argc == 2) {
        using json = nlohmann::json;
        std::ifstream file(argv[1]);
        json jsonData;
        file >> jsonData;
        size_t i = 0;
        for (const auto& item : jsonData.items()) {
            std::cout << "i " << i << std::endl;
            i++;
            input_polygons.clear();
            const auto& polygonsData = item.value();
            for (const auto& polygonData : polygonsData) {
                input_polygons.emplace_back();
                for (const auto& pointData : polygonData) {
                    double x = pointData[0];
                    double y = pointData[1];
                    input_polygons.back().push_back(CGAL::Point_2<K>(x, y));
                }
            }
            arrange_cgal_polygons(input_polygons, output);
            break;
        }
        return 0;
    } else {
        Polygon_2 rect1 = create_rectangle<typename Polygon_2::FT>(0, 0, 2, 1);
        Polygon_2 rect2 = create_rectangle<typename Polygon_2::FT>(2.2, 0, 4, 1.1);
        Polygon_2 rect3 = create_rectangle<typename Polygon_2::FT>(0, 1.2, 2, 4);
        Polygon_2 rect4 = create_rectangle<typename Polygon_2::FT>(2.2, 1.2, 6, 4);
        Polygon_2 rect5 = create_rectangle<typename Polygon_2::FT>(4.2, 0, 6, 1.1);

        input_polygons = { rect1, rect2, rect3, rect4, rect5 };
    }
    arrange_cgal_polygons(input_polygons, output);

    return 0;
}

#endif
