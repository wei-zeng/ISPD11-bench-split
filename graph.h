/*
 * a_star.h
 *
 *  Created on: Aug 16, 2019
 *      Author: wzeng
 */

#ifndef A_STAR_H_
#define A_STAR_H_

#include "baseDB.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/grid_graph.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <vector>
//#include "grDB.h"

struct VertexProp {
	// unsigned int vertexID;
	int x = -1, y = -1, z = -1;
	int real_x = -1, real_y = -1; // used in DR version where pins may be off grid
	int realPinId = -1; // for graph splitting, used in outputCSV.cpp
	bool isPin = false; // for graph splitting, used in outputCSV.cpp
	bool isRoot = false; // for graph splitting, used in outputCSV.cpp
};

struct EdgeProp {
	bool isVpin = false; // for graph splitting, used in outputCSV.cpp
	double weight = 1.0;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		VertexProp, EdgeProp> Graph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
		VertexProp, EdgeProp> DiGraph;
typedef boost::graph_traits<Graph> Traits;
typedef Graph::vertex_descriptor Vertex;
typedef Graph::vertex_iterator VertexIter;
typedef Graph::vertices_size_type VertexSize;
typedef Graph::edge_descriptor EdgeType;
typedef Graph::edge_iterator EdgeIter;
typedef Graph::edges_size_type EdgeSize;

struct vertex_hash: std::unary_function<Vertex, std::size_t> {
	std::size_t operator()(Vertex const &u) const {
		std::size_t seed = 0;
		boost::hash_combine(seed, u);
		return seed;
	}
};

struct edge_hash: std::unary_function<EdgeType, std::size_t> {
	std::size_t operator()(EdgeType const &u) const {
		std::size_t seed = 0;
		boost::hash_combine(seed, u);
		return seed;
	}
};

#endif /* A_STAR_H_ */
