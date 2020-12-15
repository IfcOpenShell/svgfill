/****************************************************************************
 * SVG fill									        					    *
 * 																			*
 * Copyright(C) 2020 AECgeeks and Bimforce								    *
 * 																		    *
 * This program is free software; you can redistribute it and/or		    *
 * modify it under the terms of the GNU Lesser General Public			    *
 * License as published by the Free Software Foundation; either			    *
 * version 3 of the License, or (at your option) any later version.		    *
 * 																		    *
 * This program is distributed in the hope that it will be useful,		    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of		    *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU	    *
 * Lesser General Public License for more details.						    *
 * 																		    *
 * You should have received a copy of the GNU Lesser General Public License *
 * along with this program; if not, write to the Free Software Foundation,  *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.	    *
 ****************************************************************************/

#include "svgfill.h"

#include <libxml/parser.h>

#include <svgpp/svgpp.hpp>
#include <svgpp/policy/xml/libxml2.hpp>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polygon_2.h>

// @todo decide on Kernel
typedef CGAL::Quotient<CGAL::MP_Float> Number_type;
typedef CGAL::Cartesian<Number_type> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::X_monotone_curve_2 Segment_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using namespace svgpp;

class Context
{
private:
	size_t depth_ = 0;
	int enabled_at_ = -1;
	svgfill::point_2 start_, xy_;

public:
	boost::optional<std::string> class_name;
	std::vector<std::vector<svgfill::line_segment_2>> segments;

	void on_enter_element(tag::element::any)
	{
		++depth_;
	}

	void on_enter_element(tag::element::g)
	{
		++depth_;
		if (enabled_at_ == -1 && !class_name.is_initialized()) {
			enabled_at_ = depth_;
			segments.emplace_back();
		}
	}

	void on_exit_element()
	{
		if (depth_-- == enabled_at_) {
			enabled_at_ = -1;
		}
	}

	template<class Str>
	void set(tag::attribute::id, Str const & s) {
	}

	template<class Str>
	void set(tag::attribute::class_, Str const & s) {
		if (enabled_at_ == -1 && (class_name.is_initialized() && s == *class_name)) {
			enabled_at_ = depth_;
			segments.emplace_back();
		}
	}

	void transform_matrix(const boost::array<double, 6> & matrix)
	{}

	void path_move_to(double x, double y, tag::coordinate::absolute)
	{
		start_ = xy_ = { x, y };
	}

	void path_line_to(double x, double y, tag::coordinate::absolute)
	{
		if (enabled_at_ != -1) {
			svgfill::point_2 next{ x, y };
			segments.back().push_back({ xy_, next });
			xy_ = next;
		}
	}

	void path_cubic_bezier_to(
		double x1, double y1,
		double x2, double y2,
		double x, double y,
		tag::coordinate::absolute) {}

	void path_quadratic_bezier_to(
		double x1, double y1,
		double x, double y,
		tag::coordinate::absolute) {}

	void path_elliptical_arc_to(
		double rx, double ry, double x_axis_rotation,
		bool large_arc_flag, bool sweep_flag,
		double x, double y,
		tag::coordinate::absolute) {}

	void path_close_subpath() {
		segments.back().push_back({ xy_, start_ });
	}

	void path_exit() {}
};

typedef
boost::mpl::set<
	// SVG Structural Elements
	tag::element::svg,
	tag::element::g,
	// SVG Shape Elements
	tag::element::circle,
	tag::element::ellipse,
	tag::element::line,
	tag::element::path,
	tag::element::polygon,
	tag::element::polyline,
	tag::element::rect
>::type processed_elements_t;

// This cryptic code just merges predefined sequences traits::shapes_attributes_by_element
// and traits::viewport_attributes with tag::attribute::transform and tag::attribute::xlink::href 
// attributes into single MPL sequence
typedef
boost::mpl::fold<
	boost::mpl::protect<
	traits::shapes_attributes_by_element
	>,
	boost::mpl::set<
	tag::attribute::id,
	tag::attribute::class_
	>::type,
	boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
>::type processed_attributes_t;

bool svgfill::svg_to_line_segments(const std::string & filename, const boost::optional<std::string>& class_name, std::vector<std::vector<line_segment_2>>& segments)
{
	Context context;
	context.class_name = class_name;

	xmlDoc* doc = xmlReadFile(filename.c_str(), NULL, 0);
	xmlNode* elem = xmlDocGetRootElement(doc);

	try {
		document_traversal<
			processed_elements<processed_elements_t>,
			processed_attributes<processed_attributes_t>
		>::load_document(elem, context);
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		return false;
	}

	segments = context.segments;

	return true;
}

namespace {
	Polygon_2 circ_to_poly(Arrangement_2::Ccb_halfedge_const_circulator circ)
	{
		Polygon_2 poly;
		auto curr = circ;
		do {
			poly.push_back(curr->source()->point());
		} while (++curr != circ);
		return poly;
	}
}

bool svgfill::line_segments_to_polygons(const std::vector<std::vector<line_segment_2>>& segments, std::vector<std::vector<polygon_2>>& polygons)
{
	for (auto& g : segments) {
		Arrangement_2 arr;

		for (auto& l : g) {
			Point_2 a(l[0][0], l[0][1]);
			Point_2 b(l[1][0], l[1][1]);
			auto ab = b - a;
			ab /= std::sqrt(CGAL::to_double(ab.squared_length()));
			// This appears to work better generally, slightly nudge the
			// end points to make sure segments intersect.
			a -= ab * 1.e-5;
			b += ab * 1.e-5;
			Segment_2 seg(a, b);
			CGAL::insert(arr, seg);
		}

		std::vector<CGAL::Polygon_2<Kernel>> ps;
		ps.reserve(std::distance(arr.faces_begin(), arr.faces_end()));

		for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) {
			const auto& f = *it;
			if (!f.is_unbounded()) {
				ps.push_back(circ_to_poly(f.outer_ccb()));
			}
		}

		// Sort polygons to have inner loops drawn over outer boundaries.
		std::sort(ps.begin(), ps.end(), [](const Polygon_2& a, const Polygon_2& b) {
			return a.area() > b.area();
		});

		polygons.emplace_back();
		polygons.back().reserve(ps.size());

		std::transform(ps.begin(), ps.end(), std::back_inserter(polygons.back()), [](const Polygon_2& p) {
			svgfill::polygon_2 p2;
			std::transform(p.vertices_begin(), p.vertices_end(), std::back_inserter(p2), [](const Point_2& pt) {
				return svgfill::point_2{
					CGAL::to_double(pt.cartesian(0)),
					CGAL::to_double(pt.cartesian(1)),
				};
			});
			return p2;
		});		
	}

	return true;
}
