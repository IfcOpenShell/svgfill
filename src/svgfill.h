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

#ifndef SVGFILL_H
#define SVGFILL_H

#include <boost/optional.hpp>

#include <array>
#include <vector>

namespace svgfill {
	typedef std::array<double, 2> point_2;
	typedef std::array<point_2, 2> line_segment_2;
	typedef std::vector<point_2> loop_2;
	struct polygon_2 {
		loop_2 boundary;
		std::vector<loop_2> inner_boundaries;
		point_2 point_inside;
	};

	enum solver {
		CARTESIAN_DOUBLE,
		CARTESIAN_QUOTIENT,
		FILTERED_CARTESIAN_QUOTIENT,
		EXACT_PREDICATES,
		EXACT_CONSTRUCTIONS
	};

	bool svg_to_line_segments(const std::string& data, const boost::optional<std::string>& class_name, std::vector<std::vector<line_segment_2>>& segments);
	bool line_segments_to_polygons(solver s, double eps, const std::vector<std::vector<line_segment_2>>& segments, std::vector<std::vector<polygon_2>>& polygons);
	bool line_segments_to_polygons(solver s, double eps, const std::vector<std::vector<line_segment_2>>& segments, std::vector<std::vector<polygon_2>>& polygons, std::function<void(float)>& progress);
	std::string polygons_to_svg(const std::vector<std::vector<polygon_2>>& polygons, bool random_color=false);
}

#endif
