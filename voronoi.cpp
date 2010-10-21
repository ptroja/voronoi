//============================================================================
// Name        : voronoi.cpp
// Author      : Piotr Trojanek
// Version     :
// Copyright   : GPL
// Description : Hello World in C, Ansi-style
//============================================================================


#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
#include <utility>                          // for std::pair
#include <exception>
#include <vector>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "CImg.h"
#include "PolyLine_Simplification.h"

using namespace cimg_library;

//! Data type for a Voronoi graph vertex
typedef struct _voronoi_vertex
{
	double x, y;
	bool occupied;
} voronoi_vertex_t;

//! Data type for a length of an graph edge
typedef double voronoi_edge_length_t;

typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, voronoi_vertex_t, voronoi_edge_length_t> Graph;
typedef boost::graph_traits <Graph>::vertex_iterator vertex_iter;
typedef boost::graph_traits <Graph>::edge_iterator edge_iter;
typedef std::vector<Graph::vertex_descriptor> graph_vertex_container;
typedef std::vector<Graph::edge_descriptor> graph_edge_container;
typedef std::pair<graph_vertex_container,graph_edge_container> graph_vertices_and_adges_container;

bool segment_cross_occupied_point(const Point & C, const Point & D, const CImg<bool> & image)
{
	// Check if the C-D segment does not cross an occupied point
	bool cross = false;
	for(unsigned int dim = 0; dim < C.size(); ++dim) {
		Point K = C;
		const int steps = (int) (D[dim] - C[dim]);
		Point Kdelta;
		if(steps) {
			Kdelta = (D - C)/(std::abs(steps));
		}
		for(int s = 0; s < std::abs(steps); ++s, K += Kdelta) {
			if(!image((int) K[0], (int) K[1])) {
				cross = true;
				break;
			}
		}
		if(cross) break;
	}

	return cross;
}

graph_vertices_and_adges_container add_path_point(const Point & C, Graph & g, const CImg<bool> & image)
{
	if(C[0] < 0 || C[1] >= image.width() || C[1] < 0 || C[1] >= image.height()) {
		throw std::logic_error("point coordinates out of map");
	}

	//! Nearest edge descriptor
	Graph::edge_descriptor ed_nearest;

	// Dummy initialize to prevent compiler warnings
	ed_nearest.m_source = (Graph::vertex_descriptor) -1;
	ed_nearest.m_target = (Graph::vertex_descriptor) -1;

	//! Nearest edge distance
	double e_nearest_distance = -1.0;

	for (std::pair <edge_iter, edge_iter> ep = boost::edges(g); ep.first != ep.second; ++ep.first) {
		const Graph::edge_descriptor ed = *ep.first;
//		if (!g[ed.m_source].occupied || !g[ed.m_target].occupied) {
			//! Coordinates of edge source
			const Point A(g[ed.m_source].x, g[ed.m_source].y);
			//! Coordinates of edge target
			const Point B(g[ed.m_target].x, g[ed.m_target].y);

			//! Segment between edge vertices
			const Segment S = { A, B };

			// Calculate distance to the edge
			Point D;
			const double d = dist_Point_to_Segment(C, S, D);

			// This is the nearest (or the first) edge
			if (d < e_nearest_distance || e_nearest_distance < 0) {
				// Check if D is within a map
				if(D[0] < 0 || D[0] >= image.width() || D[1] < 0 || D[1] >= image.height()) {
					continue;
				}

				// Check if the C-D segment does not cross an occupied point
				if(segment_cross_occupied_point(C, D, image)) {
					continue;
				}

				// OK, this is currently the best edge
				ed_nearest = ed;
				e_nearest_distance = d;
			}
//		}
	}

	if (e_nearest_distance >= 0) {
		//! Coordinates of edge source
		const Point A(g[ed_nearest.m_source].x, g[ed_nearest.m_source].y);
		//! Coordinates of edge target
		const Point B(g[ed_nearest.m_target].x, g[ed_nearest.m_target].y);

		//! Segment between edge vertices
		const Segment S = { A, B };

		Point D;
		const double C_AB = dist_Point_to_Segment(C, S, D);

		// Add the start point as a voronoi graph vertex
		voronoi_vertex_t v;

		v.occupied = false;
		v.x = C[0];
		v.y = C[1];

		const Graph::vertex_descriptor vd = add_vertex(v, g);

		//! Container for new vertices
		graph_vertex_container new_vertices;

		//! Container for new edges
		graph_edge_container new_edges;

		new_vertices.push_back(vd);

		if ((D == A).min()) {
			std::pair <Graph::edge_descriptor, bool> edge = boost::add_edge(vd, ed_nearest.m_source, g);

			// Calculate the length of the edge
			g[edge.first] = C_AB;

			new_edges.push_back(edge.first);
		} else if ((D == B).min()) {
			std::pair <Graph::edge_descriptor, bool> edge = boost::add_edge(vd, ed_nearest.m_target, g);

			// Calculate the length of the edge
			g[edge.first] = C_AB;

			new_edges.push_back(edge.first);
		} else {
			// add a cross-section of a point and a segment as a voronoi vertex
			voronoi_vertex_t v2;

			v2.occupied = false;
			v2.x = D[0];
			v2.y = D[1];

			const Graph::vertex_descriptor vd2 = add_vertex(v2, g);

			std::pair <Graph::edge_descriptor, bool> edge = boost::add_edge(vd, vd2, g);
			assert(edge.second);

			g[edge.first] = C_AB;

			// Connect a cross-section point with edge's vertices
			const std::pair <Graph::edge_descriptor, bool> edge2A = boost::add_edge(vd2, ed_nearest.m_source, g);
			assert(edge2A.second);

			if (g[ed_nearest.m_source].occupied) {
				printf("A is occupied\n");
				g[edge2A.first] = std::numeric_limits<boost::edge_bundle_type<Graph>::type>::max();
			} else {
				g[edge2A.first] = d(A,D);
			}

			const std::pair <Graph::edge_descriptor, bool> edge2B = boost::add_edge(vd2, ed_nearest.m_target, g);
			assert(edge2B.second);

			if (g[ed_nearest.m_target].occupied) {
				printf("B is occupied\n");
				g[edge2B.first] = std::numeric_limits<boost::edge_bundle_type<Graph>::type>::max();
			} else {
				g[edge2B.first] = d(B,D);
			}

			new_vertices.push_back(vd2);

			new_edges.push_back(edge.first);
			new_edges.push_back(edge2A.first);
			new_edges.push_back(edge2B.first);
		}

		return graph_vertices_and_adges_container(new_vertices, new_edges);
	}

	throw std::logic_error("unable to add vertex");
}

// Structure for point coordinates
typedef struct _point_xy {
	unsigned int x, y;
} point_xy_t;

int main(int argc, char *argv[])
{
	const std::string image_filename(argv[1]);
	const CImg <bool> orig_image(image_filename.c_str());

	// Build an image data filename
	std::string image_data_filename;
	image_data_filename += P_tmpdir;
	image_data_filename += "/";
	image_data_filename += "image_points.dat";

	// Build an voronoi graph filename
	std::string voronoi_data_filename;
	voronoi_data_filename += P_tmpdir;
	voronoi_data_filename += "/";
	voronoi_data_filename += "voronoi_points.dat";

	//*************************************************************************
	{
		// Container for coordinates of occupied grids
		std::vector<point_xy_t> points_xy;

		// Iterate over the map image to count the occupied points
		cimg_forXY(orig_image,x,y)
		{
			// printf("(%d %d)=%d\n", x, y, image(x,y));
			if (orig_image(x, y) == 0) {
				point_xy_t point;
				point.x = x;
				point.y = y;
				points_xy.push_back(point);
			}
		}

		// Open image vertices file
		std::ofstream file(image_data_filename.c_str());
		if (!file.good()) {
			perror("ostream()");
			return -1;
		}

		// Write the header (number of data dimensions)
		file << "2" << std::endl;

		// Write actual number of points
		file << points_xy.size() << std::endl;

		BOOST_FOREACH(const point_xy_t point, points_xy) {
			file << point.x << " " << point.y << std::endl;
		}

		// Exiting scope closes the data file
		file.close();
	}

	// Build the Q-hull command
	double MinFacetArea = 1.000001;
	if(argc > 2) {
		MinFacetArea = boost::lexical_cast<double>(argv[2]);
	}
	std::string command;
	command += "qvoronoi p Qt Fn PF";
	command += boost::lexical_cast<std::string>(MinFacetArea);
	for(int i = 3; i < argc; ++i) {
		std::string arg = boost::lexical_cast<std::string>(argv[i]);
		command += " ";
		command += arg;
	}
	command += " < ";
	command += image_data_filename;
	command += " > ";
	command += voronoi_data_filename;

	std::cerr << command << std::endl;

	// Call the Q-hull
	if (system(command.c_str())) {
		perror("call to qhull failed");
		return -1;
	}

	// Remove temporary input file
	if(remove(image_data_filename.c_str())) {
		perror("remove image data file()");
	}

	//*************************************************************************
	{
		// BGL data structure
		using namespace boost;

		// Open the output file
		std::ifstream data(voronoi_data_filename.c_str());

		if (!data.good()) {
			// TODO: throw
			return -1;
		}

		//! Number of dimensions
		unsigned int dim;

		//! Number of Voronoi vertices
		unsigned int number_of_vertices;

		// Read the data
		data >> dim;
		data >> number_of_vertices;

		Graph g(number_of_vertices);

		unsigned int num_of_occupied_vertices = 0;

		// Read the coordinates data
		for (std::pair <vertex_iter, vertex_iter> vp = boost::vertices(g); vp.first != vp.second; ++vp.first) {
			const Graph::vertex_descriptor vd = *vp.first;
			data >> g[vd].x;
			data >> g[vd].y;

			// Threat points outside the map as occupied
			if(g[vd].x < 0 || g[vd].x >= orig_image.width() || g[vd].y < 0 || g[vd].y >= orig_image.height()) {
				g[vd].occupied = true;
				continue;
			}

			// Check if the 3x3 neighborhood is occupied
			bool occupied = false;
			for (int xoff = -1; xoff <= 1 && !occupied; ++xoff) {
				for (int yoff = -1; yoff <= 1 && !occupied; ++yoff) {
					const int x = ((int) g[vd].x) + xoff;
					const int y = ((int) g[vd].y) + yoff;

					// Check if the pixel is within a map
					if (x >= 0 && x < (int) orig_image.width() && y >= 0 && y < (int) orig_image.height()) {
						occupied |= !orig_image(x, y);
					}
				}
			}
			g[vd].occupied = occupied;
			if(occupied) ++num_of_occupied_vertices;
		}
		printf("num_of_occupied_vertices = %d\n", num_of_occupied_vertices);

		// Read the number of vertices
		data >> number_of_vertices;

		for (unsigned int i = 0; i < number_of_vertices; ++i) {
			unsigned int number_of_neighbours;
			data >> number_of_neighbours;

			for (unsigned int j = 0; j < number_of_neighbours; ++j) {
				int neighbour;
				data >> neighbour;

				// If an edge to a defined vertex
				if (neighbour >= 0) {
					if(g[i].x < 0 || g[i].x >= orig_image.width() ||
						g[i].y < 0 || g[i].y >= orig_image.height() ||
						g[neighbour].x < 0 || g[neighbour].x >= orig_image.width() ||
						g[neighbour].y < 0 || g[neighbour].y >= orig_image.height())
						continue;
					if(segment_cross_occupied_point(Point(g[i].x, g[i].y), Point(g[neighbour].x, g[neighbour].y), orig_image))
						continue;
					std::pair <Graph::edge_descriptor, bool> edge = boost::add_edge(i, neighbour, g);

					// If this is a new edge
					if (edge.second) {
						if (g[edge.first.m_source].occupied || g[edge.first.m_target].occupied) {
							// Infinite distance between occupied vertices
							g[edge.first] = std::numeric_limits<edge_bundle_type<Graph>::type>::max();
						} else {
							// Calculate the length of the edge
							g[edge.first]
									= hypot(g[edge.first.m_target].x - g[edge.first.m_source].x, g[edge.first.m_target].y
											- g[edge.first.m_source].y);
						}
					}
				}
			}
		}

		// Remove a temporary output file
		if(remove(voronoi_data_filename.c_str())) {
			perror("remove voronoi data file()");
		}

		const bool red[] = { true, 0, 0 }, black[] = { 0, 0, 0 }, blue[] = { 0, 0, true }, green[] = {0, true, 0};

		// Data structure for debug image
		CImg <bool> visu(orig_image.width(), orig_image.height(), 1, 3, true);

		// Remove occupied vertices
		/*
		for (std::pair<vertex_iter, vertex_iter> vp = boost::vertices(g); vp.first != vp.second; ++vp.first) {
			const Graph::vertex_descriptor vd = *vp.first;
			if (g[vd].occupied) {
				boost::remove_vertex(vd, g);
			}
		}
		*/

		for (std::pair <vertex_iter, vertex_iter> vp = boost::vertices(g); vp.first != vp.second; ++vp.first) {
			const Graph::vertex_descriptor vd = *vp.first;
			//			printf("(%f,%f) -> %d\n", g[vd].x, g[vd].y, g[vd].occupied);
			visu.draw_point((int) g[vd].x, (int) g[vd].y, g[vd].occupied ? red : blue);
		}

		CImg <bool> voronoi_image(orig_image.width(), orig_image.height(), 1, 3, true);
		cimg_forXY(voronoi_image,x,y)
		{
			// printf("(%d %d)=%d\n", x, y, image(x,y));
			if (orig_image(x, y) == 0) {
				voronoi_image.draw_point(x, y, black);
			}
		}


		printf("visu: %d %d %d %d\n", visu.width(), visu.height(), visu.depth(), visu.spectrum());
		printf("voronoi_image: %d %d %d %d\n", voronoi_image.width(), voronoi_image.height(), voronoi_image.depth(), voronoi_image.spectrum());

		for (std::pair <edge_iter, edge_iter> ep = boost::edges(g); ep.first != ep.second; ++ep.first) {
			const Graph::edge_descriptor ed = *ep.first;
			//			printf("(%f,%f) -> %d\n", g[vd].x, g[vd].y, g[vd].occupied);
			voronoi_image.draw_line((int) g[ed.m_source].x, (int) g[ed.m_source].y, (int) g[ed.m_target].x, (int) g[ed.m_target].y,
					(g[ed.m_source].occupied || g[ed.m_target].occupied) ? blue : green);
		}

		int start_x = -1, start_y = -1, goal_x = -1, goal_y = -1;

		printf("initial graph vertices: %d\tedged: %d\n", g.m_vertices.size(), g.m_edges.size());

		CImgDisplay main_disp(voronoi_image, "Click a point"), draw_disp(visu, "Intensity profile");
		while (!main_disp.is_closed() && !draw_disp.is_closed()) {
			main_disp.wait();

			// printf("%d %d %d\n", main_disp.button, main_disp.mouse_y, main_disp.mouse_x);
			if (main_disp.button() && main_disp.mouse_y() >= 0 && main_disp.mouse_x()) {

				if (main_disp.button() == 1) {
					// Left button marks a start point
					start_x = main_disp.mouse_x();
					start_y = main_disp.mouse_y();
				} else if (main_disp.button() == 2) {
					// Left button marks a goal point
					goal_x = main_disp.mouse_x();
					goal_y = main_disp.mouse_y();
				}

				// Restore the original image
				CImg <bool> image(voronoi_image);

				if (start_x >= 0 && start_y >= 0) {
					image.draw_circle(start_x, start_y, 5, blue);
				}

				if(goal_x >= 0 && goal_y >= 0) {
					image.draw_circle(goal_x, goal_y, 5, red);
				}

//				printf("S %d %d (%d)\tG %d %d (%d)\n",
//						start_x, start_y, orig_image(start_x, start_y),
//						goal_x, goal_y, orig_image(goal_x, goal_y)
//				);

				if (start_x >= 0 && start_y >= 0 && goal_x >= 0 && goal_y >= 0 &&
					// Check if start or goal is not occupied
					orig_image(start_x, start_y) && orig_image(goal_x, goal_y)) {

					try {
						boost::timer cputimer;

						graph_vertices_and_adges_container start_ve_vec = add_path_point(Point(start_x, start_y), g, orig_image);
						const Graph::vertex_descriptor start_v = *(start_ve_vec.first.begin());

						graph_vertices_and_adges_container goal_ve_vec = add_path_point(Point(goal_x, goal_y), g, orig_image);
						const Graph::vertex_descriptor goal_v = *(goal_ve_vec.first.begin());

						std::vector <Graph::vertex_descriptor> p(num_vertices(g));

						boost::dijkstra_shortest_paths(g, start_v, predecessor_map(&p[0]).weight_map(get(edge_bundle, g)));

						// Simplify the polyline with Douglas-Peucker algorithm
						// http://softsurfer.com/Archive/algorithm_0205/algorithm_0205.htm

						// Containers with original and simplified vertices
						std::vector<Point> PolyLine;
						for(Graph::vertex_descriptor path_point = goal_v; path_point != p[path_point]; path_point = p[path_point]) {
							Point p(g[path_point].x, g[path_point].y);
							PolyLine.push_back(p);
						}

						printf("orignal %d", PolyLine.size());
						if (PolyLine.size() > 1) {
							std::vector<Point> sV;
							poly_simplify(5.0, PolyLine, sV);
							printf(" simplified %d", sV.size());

							// Draw the simplified results
							for(unsigned int j = 0; j < sV.size() - 1; ++j) {
								image.draw_line(
										(int) sV[j][0], (int) sV[j][1],
										(int) sV[j+1][0], (int) sV[j+1][1],
										red);
							}
						}
						printf("\ttime %.3f sec\n", cputimer.elapsed());

						// Draw the original results
						for(Graph::vertex_descriptor path_point = goal_v; path_point != p[path_point]; path_point = p[path_point]) {
							image.draw_line(
									(int) g[path_point].x, (int) g[path_point].y,
									(int) g[p[path_point]].x, (int) g[p[path_point]].y,
									black);
//							printf("(%.2f,%.2f) %d\n", g[path_point].x, g[path_point].y, g[path_point].occupied);
						}

						// Remove edges added to the graph
						BOOST_FOREACH(const Graph::edge_descriptor & ed, goal_ve_vec.second) {
							boost::remove_edge(ed, g);
						}
						BOOST_FOREACH(const Graph::edge_descriptor & ed, start_ve_vec.second) {
							boost::remove_edge(ed, g);
						}

						// Remove vertices added to the graph
						BOOST_FOREACH(const Graph::vertex_descriptor & vd, goal_ve_vec.first) {
							boost::remove_vertex(vd, g);
						}
						BOOST_FOREACH(const Graph::vertex_descriptor & vd, start_ve_vec.first) {
							boost::remove_vertex(vd, g);
						}

					} catch (std::logic_error & e) {
						std::cerr << "error: " << e.what() << std::endl;
					}
				}

				// Display results
				main_disp.display(image);
			}
		}
	}

	return 0;
}
