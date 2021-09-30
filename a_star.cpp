#include "graph.h"
#include "grDB.h"
#include <map>
#include <vector>
#include <set>
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <chrono>
#include <cmath>

// Manhattan distance heuristic
template<typename Graph, typename CostType>
class distance_heuristic: public boost::astar_heuristic<Graph, CostType> {
public:
	distance_heuristic(const unordered_set<Vertex> &goals, Graph graph) :
			m_graph(graph), m_goals(goals) {
	}
	CostType operator()(Vertex u) {
		CostType result = INFINITY;
		for (Vertex m_goal : m_goals) {
			CostType dx = m_graph[m_goal].x - m_graph[u].x;
			CostType dy = m_graph[m_goal].y - m_graph[u].y;
			CostType dz = m_graph[m_goal].z - m_graph[u].z;
			CostType tmp = fabs(dx) + fabs(dy) + fabs(dz);
			if (tmp < result)
				result = tmp;
		}
		return result;
	}
private:
	Graph m_graph;
	unordered_set<Vertex> m_goals;
};

struct found_goal {
	Vertex goal;
};

struct astar_goal_visitor: public boost::default_astar_visitor {
	astar_goal_visitor(const unordered_set<Vertex> &goals) :
			m_goals(goals) {
	}
	void examine_vertex(Vertex u, const Graph &g) {
		if (m_goals.find(u) != m_goals.end()) {
			throw found_goal( { u });
		}
	}
private:
	unordered_set<Vertex> m_goals;
	// boost::unordered_map<Vertex, Vertex, vertex_hash> m_predecessor;
};

std::pair<double, vector<Vertex>> RoutingDB::a_star(Graph &graph3D,
		const Vertex &s, const unordered_set<Vertex> &goals) const {
	using namespace boost;
	typedef boost::unordered_map<Vertex, Vertex, vertex_hash> pred_map;
	pred_map predecessor;
	boost::associative_property_map<pred_map> pred_pmap(predecessor);

	typedef boost::unordered_map<Vertex, double, vertex_hash> dist_map;
	dist_map distance;
	boost::associative_property_map<dist_map> dist_pmap(distance);

	typedef boost::unordered_map<Vertex, double, vertex_hash> cost_map;
	cost_map cost;
	boost::associative_property_map<cost_map> cost_pmap(cost);

	typedef boost::unordered_map<EdgeType, double, edge_hash> wt_map;
	wt_map weight;
	EdgeIter e, eend;
	for (tie(e, eend) = edges(graph3D); e != eend; ++e) {
		weight[*e] = graph3D[*e].weight;
	}
	Graph::out_edge_iterator oe, oeend;
	tie(oe, oeend) = out_edges(s, graph3D);
	boost::associative_property_map<wt_map> wt_pmap(weight);

	astar_goal_visitor visitor(goals);

	auto heur = distance_heuristic<Graph, double>(goals, graph3D);
	vector<Vertex> vecPath;
	try {
		VertexIter ui, ui_end;
		for (boost::tie(ui, ui_end) = vertices(graph3D); ui != ui_end; ++ui) {
			put(dist_pmap, *ui, INFINITY);
			put(cost_pmap, *ui, INFINITY);
			put(pred_pmap, *ui, *ui);
			visitor.initialize_vertex(*ui, graph3D);
		}
		put(dist_pmap, s, 0);
		put(cost_pmap, s, heur(s));
		boost::astar_search_no_init_tree(graph3D, s, heur,
				weight_map(wt_pmap).predecessor_map(pred_pmap).distance_map(
						dist_pmap).visitor(visitor));
		// printf("No feasible path.\n");
		return make_pair(INFINITY, vecPath);
	} catch (found_goal &fg) {
		// Walk backwards from the goal through the predecessor chain adding
		// vertices to the solution path.
		//printf("Path found: ");
		//isFeasible = true;
		double wl = 0;
		for (Vertex u = fg.goal; u != s; u = predecessor[u]) {
			//	printf("(%d,%d,%d)-", graph3D[u].x, graph3D[u].y, graph3D[u].z);
			// auto r = edge(u, predecessor[u], graph3D);
			// assert(r.second);
			vecPath.push_back(u);
			//	graph3D[r.first].isPath = true;
			++wl;
		}
		vecPath.push_back(s);
		// printf("(%d,%d,%d)\n", graph3D[s].x, graph3D[s].y, graph3D[s].z);
		// cout << "WL = " << wl << endl;
		return make_pair(wl, vecPath);
	}
}

unordered_set<Vertex> RoutingDB::build_goals(const Vpin &vp, const int xlLim,
		const int xuLim, const int ylLim, const int yuLim,
		const int zlLim) const {
	using namespace boost;
	std::unordered_set<Vertex> goals;
	int X = xuLim - xlLim + 1, Y = yuLim - ylLim + 1;
	for (auto gv : vp.gnetGrids) {
		goals.insert(
				(gv._z - zlLim) * Y * X + (gv._y - ylLim) * X
						+ (gv._x - xlLim));
	}
	return goals;
}

Graph RoutingDB::build_graph(const size_t gnetId,
		const std::unordered_set<Gcell, HashGcell3d> &bannedPts,
		const int xlLim, const int xuLim, const int ylLim, const int yuLim,
		const int zlLim, const int zuLim) const {
	using namespace boost;
	Graph graph3D;
	// vector<double> capacities;
	// vector<Location> locations;
	// boost::property_map<Graph, boost::edge_index_t>::type edge_index_map = get(boost::edge_index, graph2D);
	// boost::property_map<Graph, boost::vertex_index_t>::type vertex_index_map = get(boost::vertex_index, graph2D);
	int X = xuLim - xlLim + 1, Y = yuLim - ylLim + 1;
	for (int z = zlLim; z <= zuLim; ++z) {
		for (int y = ylLim; y <= yuLim; ++y) {
			for (int x = xlLim; x <= xuLim; ++x) {
				Vertex v = add_vertex(graph3D);
				graph3D[v].x = x;
				graph3D[v].y = y;
				graph3D[v].z = z;
			}
		}
	}
	int cap, dem, blk;
	bool vacant, banned;
	for (int x = xlLim; x <= xuLim; ++x) {
		for (int y = ylLim; y <= yuLim; ++y) {
			if (x > xlLim) {
				for (int z = zlLim; z <= zuLim; ++z) {
					if (getRoutingDir(z) == H) {
						size_t edgeId = findEdge(x - 1, y, z);
						cap = getEdgeCap(edgeId);
						dem = getEdgeDemand(edgeId);
						blk = getEdgeBlk(edgeId);
						vacant = _gnets[gnetId]._occupiedEdges.find(edgeId)
								== _gnets[gnetId]._occupiedEdges.end();
						banned = (bannedPts.find(Gcell(x - 1, y, z))
								!= bannedPts.end()
								|| bannedPts.find(Gcell(x, y, z))
										!= bannedPts.end());
					} else {
						cap = 0;
						dem = 0;
						blk = 0;
						vacant = false;
						banned = false;
					}
					if (cap > dem + blk && vacant && !banned) {
						auto r = add_edge(
								vertex(
										(z - zlLim) * Y * X + (y - ylLim) * X
												+ (x - xlLim), graph3D),
								vertex(
										(z - zlLim) * Y * X + (y - ylLim) * X
												+ (x - 1 - xlLim), graph3D),
								graph3D);
						assert(r.second);
					}
				}
			}
			if (y > ylLim) {
				for (int z = zlLim; z <= zuLim; ++z) {
					if (getRoutingDir(z) == V) {
						size_t edgeId = findEdge(x, y - 1, z);
						cap = getEdgeCap(edgeId);
						dem = getEdgeDemand(edgeId);
						blk = getEdgeBlk(edgeId);
						vacant = _gnets[gnetId]._occupiedEdges.find(edgeId)
								== _gnets[gnetId]._occupiedEdges.end();
						banned = (bannedPts.find(Gcell(x, y - 1, z))
								!= bannedPts.end()
								|| bannedPts.find(Gcell(x, y, z))
										!= bannedPts.end());
					} else {
						cap = 0;
						dem = 0;
						blk = 0;
						vacant = false;
						banned = false;
					}
					if (cap > dem + blk && vacant && !banned) {
						auto r = add_edge(
								vertex(
										(z - zlLim) * Y * X + (y - ylLim) * X
												+ (x - xlLim), graph3D),
								vertex(
										(z - zlLim) * Y * X
												+ (y - 1 - ylLim) * X
												+ (x - xlLim), graph3D),
								graph3D);
						assert(r.second);
					}
				}
			}
			for (int z = zlLim + 1; z <= zuLim; ++z) {
				size_t viaId = findVia(x, y, z - 1);
				cap = getViaCap(viaId);
				dem = getViaDemand(viaId);
				vacant = _gnets[gnetId]._occupiedVias.find(viaId)
						== _gnets[gnetId]._occupiedVias.end();
				banned = (bannedPts.find(Gcell(x, y, z - 1)) != bannedPts.end()
						|| bannedPts.find(Gcell(x, y, z)) != bannedPts.end());
				if (cap > dem && vacant && !banned) {
					auto r = add_edge(
							vertex(
									(z - zlLim) * Y * X + (y - ylLim) * X
											+ (x - xlLim), graph3D),
							vertex(
									(z - 1 - zlLim) * Y * X + (y - ylLim) * X
											+ (x - xlLim), graph3D), graph3D);
					assert(r.second);
				}
			}
		}
	}
	return graph3D;
}

bool comp_x(Gcell i, Gcell j) {
	return i._x < j._x;
}
bool comp_y(Gcell i, Gcell j) {
	return i._y < j._y;
}

std::tuple<double, vector<int>, vector<int>, int,
		std::unordered_set<Gcell, HashGcell3d>> RoutingDB::rerouteNet(
		const Layout &layout, Vpin &vp,
		const std::unordered_set<Gcell, HashGcell3d> &bannedPts,
		const int zlLim, const int zuLim, const int maxMargin,
		const bool singleMargin, const bool writeToDB, const bool updateRoot) {
	// Connect one v-pin to its remaining net
	using namespace boost;
	auto terminals = vp.gnetVertices;
	std::unordered_set<Gcell, HashGcell3d> addlBannedPts;
	terminals.insert(Gcell(vp.xCoord, vp.yCoord, vp.zCoord));
	int margin = singleMargin ? maxMargin : 0;
	int next_margin = 0;
	double wl = INFINITY;
	bool isFeasible = false;
	double minmhd = INFINITY;
	std::pair<double, vector<Vertex>> res;
	Graph g;
	std::vector<int> addedEdges, addedVias;
	int viaID = findVia(vp.xCoord, vp.yCoord, vp.zCoord);
	if (getViaDemand(viaID) >= getViaCap(viaID)) {
		return make_tuple(wl, addedEdges, addedVias, margin, addlBannedPts); // infeasible
	}
	// add new vpin
	// addViaDemand(viaID);
	// addedVias.push_back(viaID);
	while (margin <= maxMargin) {
// cout << "margin = " << margin << endl;
		int xlLim = max(
				min_element(terminals.begin(), terminals.end(), comp_x)->_x
						- margin, 0);
		int xuLim = min(
				max_element(terminals.begin(), terminals.end(), comp_x)->_x
						+ margin, (int) layout._numTilesX - 1);
		int ylLim = max(
				min_element(terminals.begin(), terminals.end(), comp_y)->_y
						- margin, 0);
		int yuLim = min(
				max_element(terminals.begin(), terminals.end(), comp_y)->_y
						+ margin, (int) layout._numTilesY - 1);
		g = build_graph(vp.gnetID, bannedPts, xlLim, xuLim, ylLim, yuLim, zlLim,
				zuLim);
		/*for (auto v : g.m_vertices) {
		 cerr << "Vertex: " << v.m_property.x << " "
		 << v.m_property.y << " " << v.m_property.z << endl;
		 }
		 for (auto e : g.m_edges) {
		 cerr << "Edge: " << g[e.m_source].id << " " << g[e.m_target].id
		 << endl;
		 }*/

		int X = xuLim - xlLim + 1, Y = yuLim - ylLim + 1;
		Vertex start = (vp.zCoord - zlLim) * Y * X + (vp.yCoord - ylLim) * X
				+ (vp.xCoord - xlLim);
		std::unordered_set<Vertex> goals = build_goals(vp, xlLim, xuLim, ylLim,
				yuLim, zlLim);
		if (minmhd == INFINITY) {
			for (auto goal : goals) {
				double mhd = fabs(g[start].x - g[goal].x)
						+ fabs(g[start].y - g[goal].y)
						+ fabs(g[start].z - g[goal].z);
				if (mhd < minmhd)
					minmhd = mhd;
			}
		}
		res = a_star(g, start, goals);
		if (res.first < INFINITY) { // found a-star path
			wl = res.first;
// cout << "wl = " << wl << endl;
			isFeasible = true;
			next_margin = std::max(next_margin,
					std::min(maxMargin, (int) ((res.first - minmhd) / 2.0)));
		} else { // not found
			next_margin = std::max(next_margin,
					std::min(maxMargin, margin * 2 + 1));
// cout << "next_margin = " << next_margin << endl;
		}
		if (next_margin > margin)
			margin = next_margin;
		else
			break;
	}
	if (!isFeasible)
		return make_tuple(INFINITY, addedEdges, addedVias, margin,
				addlBannedPts);
	// fill addedEdges / addedVias
	auto vecPath = res.second;
	if (writeToDB) {
		vp.wlToL1 += vecPath.size();
	}
	if (vecPath.size() > 0)
		addlBannedPts.insert(
				Gcell(g[vecPath[0]].x, g[vecPath[0]].y, g[vecPath[0]].z));
	for (int i = 1; i < static_cast<int>(vecPath.size()); i++) {
		addlBannedPts.insert(
				Gcell(g[vecPath[i]].x, g[vecPath[i]].y, g[vecPath[i]].z));
		if (g[vecPath[i]].z == g[vecPath[i - 1]].z) {
			if (_dirLayers[g[vecPath[i]].z] == H) {
				assert(g[vecPath[i]].y == g[vecPath[i - 1]].y);
				for (int x = min(g[vecPath[i]].x, g[vecPath[i - 1]].x);
						x < max(g[vecPath[i]].x, g[vecPath[i - 1]].x); x++) {
					if (writeToDB) {
						int edgeId = findEdge(x, g[vecPath[i]].y,
								g[vecPath[i]].z);
						if (!(_gnets[vp.gnetID].findCleanEdge(edgeId))) {
							_gnets[vp.gnetID].setCleanEdge(edgeId,
									_edges[edgeId]);
						}
						if (!(_gnets[vp.gnetID].findRUEdge(edgeId))) {
							_gnets[vp.gnetID].setRUEdge(edgeId, _edges[edgeId]);
						}
						restoreEdge(edgeId);
						_gnets[vp.gnetID]._occupiedEdges.insert(edgeId);
						_dirty = true;
					}
					addedEdges.push_back(
							findEdge(x, g[vecPath[i]].y, g[vecPath[i]].z));
				}
			} else if (_dirLayers[g[vecPath[i]].z] == V) {
				assert(g[vecPath[i]].x == g[vecPath[i - 1]].x);
				for (int y = min(g[vecPath[i]].y, g[vecPath[i - 1]].y);
						y < max(g[vecPath[i]].y, g[vecPath[i - 1]].y); y++) {
					if (writeToDB) {
						int edgeId = findEdge(g[vecPath[i]].x, y,
								g[vecPath[i]].z);
						if (!(_gnets[vp.gnetID].findCleanEdge(edgeId))) {
							_gnets[vp.gnetID].setCleanEdge(edgeId,
									_edges[edgeId]);
						}
						if (!(_gnets[vp.gnetID].findRUEdge(edgeId))) {
							_gnets[vp.gnetID].setRUEdge(edgeId, _edges[edgeId]);
						}
						restoreEdge(edgeId);
						_gnets[vp.gnetID]._occupiedEdges.insert(edgeId);
						_dirty = true;
					}
					addedEdges.push_back(
							findEdge(g[vecPath[i]].x, y, g[vecPath[i]].z));
				}
			} else
				assert(0);
		} else {
			assert(g[vecPath[i]].y == g[vecPath[i - 1]].y);
			assert(g[vecPath[i]].x == g[vecPath[i - 1]].x);
			for (int z = min(g[vecPath[i]].z, g[vecPath[i - 1]].z);
					z < max(g[vecPath[i]].z, g[vecPath[i - 1]].z); z++) {
				if (writeToDB) {
					int viaId = findVia(g[vecPath[i]].x, g[vecPath[i]].y, z);
					if (!(_gnets[vp.gnetID].findCleanVia(viaId))) {
						_gnets[vp.gnetID].setCleanVia(viaId, _vias[viaId]);
					}
					if (!(_gnets[vp.gnetID].findRUVia(viaId))) {
						_gnets[vp.gnetID].setRUVia(viaId, _vias[viaId]);
					}
					restoreVia(viaId);
					_gnets[vp.gnetID]._occupiedVias.insert(viaId);
					_dirty = true;
				}
				addedVias.push_back(
						findVia(g[vecPath[i]].x, g[vecPath[i]].y, z));
			}
		}
	}
	if (writeToDB) {
		auto &vecPath = res.second;
		// Erase points in the middle in the same direction along the path
		vector<bool> validPt;
		validPt.assign(vecPath.size(), true);
		char prevDir = 0; // 0 = default, 1 = +x, 2 = -x, 3 = +y, 4 = -y, 5 = +z, 6 = -z, 7 = other
		if (vecPath.size() > 2) {
			for (size_t i = 1; i < vecPath.size(); ++i) {
				if (g[vecPath[i]].x > g[vecPath[i - 1]].x
						&& g[vecPath[i]].y == g[vecPath[i - 1]].y
						&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
					if (prevDir == 1) {
						validPt[i - 1] = false;
					} else {
						prevDir = 1;
					}
				} else if (g[vecPath[i]].x < g[vecPath[i - 1]].x
						&& g[vecPath[i]].y == g[vecPath[i - 1]].y
						&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
					if (prevDir == 2) {
						validPt[i - 1] = false;
					} else {
						prevDir = 2;
					}
				} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
						&& g[vecPath[i]].y > g[vecPath[i - 1]].y
						&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
					if (prevDir == 3) {
						validPt[i - 1] = false;
					} else {
						prevDir = 3;
					}
				} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
						&& g[vecPath[i]].y < g[vecPath[i - 1]].y
						&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
					if (prevDir == 4) {
						validPt[i - 1] = false;
					} else {
						prevDir = 4;
					}
				} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
						&& g[vecPath[i]].y == g[vecPath[i - 1]].y
						&& g[vecPath[i]].z > g[vecPath[i - 1]].z) {
					if (prevDir == 5) {
						// validPt[i - 1] = false;
					} else {
						prevDir = 5;
					}
				} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
						&& g[vecPath[i]].y == g[vecPath[i - 1]].y
						&& g[vecPath[i]].z < g[vecPath[i - 1]].z) {
					if (prevDir == 6) {
						// validPt[i - 1] = false;
					} else {
						prevDir = 6;
					}
				} else {
					prevDir = 7;
				}
			}
			vector<Vertex> vecPath_tmp;
			for (size_t i = 0; i < vecPath.size(); ++i) {
				if (validPt[i]) {
					vecPath_tmp.push_back(vecPath[i]);
				}
			}
			vecPath.swap(vecPath_tmp);
		}
		Gnet &gnet = _gnets[vp.gnetID];
		gnet._gWires.reserve(gnet._gWires.size() + vecPath.size() + 1);
		int addedVertex = -1; // gWireID of the added vertex, -1 if n/a.

		vector<GlobalWire>::iterator gWireItr;
// check if there is a removed vertex due to broken wires and if so, add it back
		if (updateRoot) { // this can happen only when updateRoot == true
			vector<vector<GlobalWire>::iterator> rootItrs;
			for (gWireItr = gnet._gWires.begin(); gWireItr < gnet._gWires.end();
					++gWireItr) {
				if (gWireItr->_pWireId == -1
						&& gWireItr != gnet._gWires.begin() + gnet._root) {
					rootItrs.push_back(gWireItr);
				}
			}
			auto s = vp.gnetVertices;
			for (auto gw : gnet._gWires) {
//	cout << gw._x << " " << gw._y << " " << gw._z << endl;
				s.erase(Gcell(gw._x, gw._y, gw._z));
			}
//		cout << "s.size() = " << s.size() << endl;
			if (s.size() > 1) {
				cout << s.size() << endl;
				exit(-2);
			}
			if (s.size() == 1) { // have to make a new gWire that ALL roots (except gWire[0]) point to.
				gnet._gWires.push_back(
						GlobalWire(s.begin()->_x, s.begin()->_y, s.begin()->_z,
								-1, -1));
				addedVertex = gnet._gWires.size() - 1;
				for (auto &rootItr : rootItrs) {
					rootItr->_pWireId = addedVertex;
				}
			}
		}
// Identify goal type
		for (gWireItr = gnet._gWires.begin(); gWireItr < gnet._gWires.end();
				++gWireItr) {
			if (gWireItr->_x == g[vecPath[0]].x
					&& gWireItr->_y == g[vecPath[0]].y
					&& gWireItr->_z == g[vecPath[0]].z) {
				break; // Case 0: goal is one of endpoints of a gWire
			}
		}

		if (gWireItr == gnet._gWires.end()) { // Case 1: goal in the middle of a gWire
			for (gWireItr = gnet._gWires.begin(); gWireItr < gnet._gWires.end();
					++gWireItr) {
				if (gWireItr->_pWireId != -1) {
					int &x1 = gWireItr->_x, &y1 = gWireItr->_y, &z1 =
							gWireItr->_z;
					int &x2 = gnet._gWires[gWireItr->_pWireId]._x, &y2 =
							gnet._gWires[gWireItr->_pWireId]._y, &z2 =
							gnet._gWires[gWireItr->_pWireId]._z;
					int &x0 = g[vecPath[0]].x, &y0 = g[vecPath[0]].y, &z0 =
							g[vecPath[0]].z;
					if ((x1 != x2 && y0 == y1 && z0 == z1 && min(x1, x2) < x0
							&& x0 < max(x1, x2))
							|| (y1 != y2 && x0 == x1 && z0 == z1
									&& min(y1, y2) < y0 && y0 < max(y1, y2))
							|| (z1 != z2 && x0 == x1 && y0 == y1
									&& min(z1, z2) < z0 && z0 < max(z1, z2))) {
						// make a new gWire in the middle
						gnet._gWires.push_back(
								GlobalWire(x0, y0, z0, gWireItr->_pWireId, -1));
						gWireItr->_pWireId = gnet._gWires.size() - 1;
						gWireItr = gnet._gWires.end() - 1;
						break;
					}
				}
			}
		}
		assert(gWireItr < gnet._gWires.end());
		if (updateRoot) {
// Mode 2: move the root, make the root point to the new vertex, update the root
// search for all roots due to broken wires
			auto rootItr = gWireItr;
			vector<int> gWireIDsToReverse;
			if (rootItr < gnet._gWires.end()) { // Case 0 or Case 1
				while (rootItr->_pWireId != -1) {
					gWireIDsToReverse.push_back(rootItr - gnet._gWires.begin());
					rootItr = gnet._gWires.begin() + rootItr->_pWireId;
				}
				gWireIDsToReverse.push_back(rootItr - gnet._gWires.begin());
			}

			gnet._gWires[gWireIDsToReverse[0]]._pWireId = -1;
			for (auto itr = gWireIDsToReverse.begin() + 1;
					itr < gWireIDsToReverse.end(); ++itr) {
// Reverse the directions of each gWire
				gnet._gWires[(*itr)]._pWireId = *(itr - 1);
			}
			if (gnet._gWires[gnet._root]._pWireId != -1) {
				gnet._root = gWireIDsToReverse[0];
				cout << "<I> Root changed to: " << gnet._root << endl;
			}
// let the new gWire point to the last gwire in the list to reverse,
// finishing the reversing process.
//if (addedVertex != -1) {
//	gnet._gWires[addedVertex]._pWireId =
//			*(gWireIDsToReverse.end() - 1);
//}

// expose the new root and prepare for adding the new path.
			gWireItr = gnet._gWires.begin() + gWireIDsToReverse[0];

			for (auto itr = vecPath.begin() + 1; itr < vecPath.end(); ++itr) {
// Add edge (*itr, *(itr-1)) to DB (gnet/addEdgeDemand)
				gnet._gWires.push_back(
						GlobalWire(g[*itr].x, g[*itr].y, g[*itr].z, -1, -1));
				gWireItr->_pWireId = gnet._gWires.size() - 1;
				gWireItr = gnet._gWires.end() - 1;
			}
		} else {
// Mode 1: make the new vertex point to the root, keep the root
			int pWireId = gWireItr - gnet._gWires.begin();
			for (auto itr = vecPath.begin() + 1; itr < vecPath.end(); ++itr) {
// Add edge (*itr, *(itr+1)) to DB (gnet/addEdgeDemand)
				gnet._gWires.push_back(
						GlobalWire(g[*itr].x, g[*itr].y, g[*itr].z, pWireId,
								-1));
				pWireId = gnet._gWires.size() - 1;
			}
		}
	}
	return make_tuple(wl, addedEdges, addedVias, margin, addlBannedPts);
}

std::tuple<double, vector<int>, vector<int>, int> RoutingDB::rerouteNet(
		const Layout &layout, const Vpin &vp1, const Vpin &vp2, const int zlLim,
		const int zuLim, const int maxMargin, const bool singleMargin,
		const bool writeToDB) {
// Connect two v-pins
	using namespace boost;
	int margin = singleMargin ? maxMargin : 0;
	int next_margin = 0;
	double wl = INFINITY;
	bool isFeasible = false;
	double minmhd = INFINITY;
	std::pair<double, vector<Vertex>> res;
	std::vector<int> addedEdges, addedVias;
	Graph g;
	std::unordered_set<Gcell, HashGcell3d> bannedPts; // empty set
	while (margin <= maxMargin) {
// cout << "[2vp] margin = " << margin << endl;
		int xlLim = max(min(vp1.xCoord, vp2.xCoord) - margin, 0);
		int xuLim = min(max(vp1.xCoord, vp2.xCoord) + margin,
				(int) layout._numTilesX - 1);
		int ylLim = max(min(vp1.yCoord, vp2.yCoord) - margin, 0);
		int yuLim = min(max(vp1.yCoord, vp2.yCoord) + margin,
				(int) layout._numTilesY - 1);
		g = build_graph(vp1.gnetID, bannedPts, xlLim, xuLim, ylLim, yuLim,
				zlLim, zuLim);
		int X = xuLim - xlLim + 1, Y = yuLim - ylLim + 1;
		Vertex start = (vp1.zCoord + 1 - zlLim) * Y * X
				+ (vp1.yCoord - ylLim) * X + (vp1.xCoord - xlLim);
		Vertex goal = (vp2.zCoord + 1 - zlLim) * Y * X
				+ (vp2.yCoord - ylLim) * X + (vp2.xCoord - xlLim);
		if (minmhd == INFINITY) {
			minmhd = fabs(g[start].x - g[goal].x) + fabs(g[start].y - g[goal].y)
					+ fabs(g[start].z - g[goal].z);
		}
		std::unordered_set<Vertex> goals = { goal };
		res = a_star(g, start, goals);
		if (res.first < INFINITY) { // found a-star path
			wl = res.first;
// cout << "wl = " << wl << endl;
			isFeasible = true;
			next_margin = std::max(next_margin,
					std::min(maxMargin, (int) ((res.first - minmhd) / 2.0)));
		} else { // not found
			next_margin = std::max(next_margin,
					std::min(maxMargin, margin * 2 + 1));
// cout << "next_margin = " << next_margin << endl;
		}
		if (next_margin > margin)
			margin = next_margin;
		else
			break;
	}
	if (!isFeasible)
		return make_tuple(INFINITY, addedEdges, addedVias, margin);
// add the new vpin
	if (writeToDB) {
		int viaId = findVia(vp1.xCoord, vp1.yCoord, vp1.zCoord);
		if (!(_gnets[vp1.gnetID].findCleanVia(viaId))) {
			_gnets[vp1.gnetID].setCleanVia(viaId, _vias[viaId]);
		}
		if (!(_gnets[vp1.gnetID].findRUVia(viaId))) {
			_gnets[vp1.gnetID].setRUVia(viaId, _vias[viaId]);
		}
		restoreVia(viaId);
		_dirty = true;
	}
	addedVias.push_back(findVia(vp1.xCoord, vp1.yCoord, vp1.zCoord));

	auto &vecPath = res.second;
	// Erase points in the middle in the same direction along the path
	vector<bool> validPt;
	validPt.assign(vecPath.size(), true);
	char prevDir = 0; // 0 = default, 1 = +x, 2 = -x, 3 = +y, 4 = -y, 5 = +z, 6 = -z, 7 = other
	if (vecPath.size() > 2) {
		for (size_t i = 1; i < vecPath.size(); ++i) {
			if (g[vecPath[i]].x > g[vecPath[i - 1]].x
					&& g[vecPath[i]].y == g[vecPath[i - 1]].y
					&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
				if (prevDir == 1) {
					validPt[i - 1] = false;
				} else {
					prevDir = 1;
				}
			} else if (g[vecPath[i]].x < g[vecPath[i - 1]].x
					&& g[vecPath[i]].y == g[vecPath[i - 1]].y
					&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
				if (prevDir == 2) {
					validPt[i - 1] = false;
				} else {
					prevDir = 2;
				}
			} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
					&& g[vecPath[i]].y > g[vecPath[i - 1]].y
					&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
				if (prevDir == 3) {
					validPt[i - 1] = false;
				} else {
					prevDir = 3;
				}
			} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
					&& g[vecPath[i]].y < g[vecPath[i - 1]].y
					&& g[vecPath[i]].z == g[vecPath[i - 1]].z) {
				if (prevDir == 4) {
					validPt[i - 1] = false;
				} else {
					prevDir = 4;
				}
			} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
					&& g[vecPath[i]].y == g[vecPath[i - 1]].y
					&& g[vecPath[i]].z > g[vecPath[i - 1]].z) {
				if (prevDir == 5) {
					// validPt[i - 1] = false;
				} else {
					prevDir = 5;
				}
			} else if (g[vecPath[i]].x == g[vecPath[i - 1]].x
					&& g[vecPath[i]].y == g[vecPath[i - 1]].y
					&& g[vecPath[i]].z < g[vecPath[i - 1]].z) {
				if (prevDir == 6) {
					// validPt[i - 1] = false;
				} else {
					prevDir = 6;
				}
			} else {
				prevDir = 7;
			}
		}
		vector<Vertex> vecPath_tmp;
		for (size_t i = 0; i < vecPath.size(); ++i) {
			if (validPt[i]) {
				vecPath_tmp.push_back(vecPath[i]);
			}
		}
		vecPath.swap(vecPath_tmp);
	}
	for (int i = 1; i < static_cast<int>(vecPath.size()); i++) {
		if (g[vecPath[i]].z == g[vecPath[i - 1]].z) {
			if (_dirLayers[g[vecPath[i]].z] == H) {
				assert(g[vecPath[i]].y == g[vecPath[i - 1]].y);
				for (int x = min(g[vecPath[i]].x, g[vecPath[i - 1]].x);
						x < max(g[vecPath[i]].x, g[vecPath[i - 1]].x); x++) {
					if (writeToDB) {
						int edgeId = findEdge(x, g[vecPath[i]].y,
								g[vecPath[i]].z);
						if (!(_gnets[vp1.gnetID].findCleanEdge(edgeId))) {
							_gnets[vp1.gnetID].setCleanEdge(edgeId,
									_edges[edgeId]);
						}
						if (!(_gnets[vp1.gnetID].findRUEdge(edgeId))) {
							_gnets[vp1.gnetID].setRUEdge(edgeId,
									_edges[edgeId]);
						}
						restoreEdge(edgeId);
						_dirty = true;
					}
					addedEdges.push_back(
							findEdge(x, g[vecPath[i]].y, g[vecPath[i]].z));

				}
			} else if (_dirLayers[g[vecPath[i]].z] == V) {
				assert(g[vecPath[i]].x == g[vecPath[i - 1]].x);
				for (int y = min(g[vecPath[i]].y, g[vecPath[i - 1]].y);
						y < max(g[vecPath[i]].y, g[vecPath[i - 1]].y); y++) {
					if (writeToDB) {
						int edgeId = findEdge(g[vecPath[i]].x, y,
								g[vecPath[i]].z);
						if (!(_gnets[vp1.gnetID].findCleanEdge(edgeId))) {
							_gnets[vp1.gnetID].setCleanEdge(edgeId,
									_edges[edgeId]);
						}
						if (!(_gnets[vp1.gnetID].findRUEdge(edgeId))) {
							_gnets[vp1.gnetID].setRUEdge(edgeId,
									_edges[edgeId]);
						}
						restoreEdge(edgeId);
						_dirty = true;
					}
					addedEdges.push_back(
							findEdge(g[vecPath[i]].x, y, g[vecPath[i]].z));
				}
			} else
				assert(0);
		} else {
			assert(g[vecPath[i]].y == g[vecPath[i - 1]].y);
			assert(g[vecPath[i]].x == g[vecPath[i - 1]].x);
			for (int z = min(g[vecPath[i]].z, g[vecPath[i - 1]].z);
					z < max(g[vecPath[i]].z, g[vecPath[i - 1]].z); z++) {
				if (writeToDB) {
					int viaId = findVia(g[vecPath[i]].x, g[vecPath[i]].y, z);
					if (!(_gnets[vp1.gnetID].findCleanVia(viaId))) {
						_gnets[vp1.gnetID].setCleanVia(viaId, _vias[viaId]);
					}
					if (!(_gnets[vp1.gnetID].findRUVia(viaId))) {
						_gnets[vp1.gnetID].setRUVia(viaId, _vias[viaId]);
					}
					restoreVia(findVia(g[vecPath[i]].x, g[vecPath[i]].y, z));
					_dirty = true;
				}
				addedVias.push_back(
						findVia(g[vecPath[i]].x, g[vecPath[i]].y, z));
			}
		}
	}
	if (writeToDB) {
		auto vecPath = res.second;
		Gnet &gnet = _gnets[vp1.gnetID];
		gnet._gWires.reserve(gnet._gWires.size() + vecPath.size());
		vector<GlobalWire>::iterator gWireItr1, gWireItr2;
		// start
		for (gWireItr1 = gnet._gWires.begin(); gWireItr1 < gnet._gWires.end();
				++gWireItr1) {
			if (gWireItr1->_x == g[vecPath.back()].x
					&& gWireItr1->_y == g[vecPath.back()].y
					&& gWireItr1->_z == g[vecPath.back()].z - 1) {
				break;
			}
		}
		// goal
		for (gWireItr2 = gnet._gWires.begin(); gWireItr2 < gnet._gWires.end();
				++gWireItr2) {
			if (gWireItr2->_x == g[vecPath[0]].x
					&& gWireItr2->_y == g[vecPath[0]].y
					&& gWireItr2->_z == g[vecPath[0]].z) {
				break;
			}
		}
		assert((gWireItr1->_pWireId == -1) != (gWireItr2->_pWireId == -1));
		if (gWireItr2->_pWireId == -1) {
// cout << "Mode 2" << endl;
// Mode 2
			if (vecPath.size() == 1) {
				gWireItr2->_pWireId = gWireItr1 - gnet._gWires.begin();
			} else {
				for (auto itr = vecPath.begin() + 1; itr < vecPath.end() - 1;
						++itr) {
// Add edge (*itr, *(itr+1)) to DB (gnet/addEdgeDemand)
					gnet._gWires.push_back(
							GlobalWire(g[*itr].x, g[*itr].y, g[*itr].z, -1,
									-1));
					gWireItr2->_pWireId = gnet._gWires.size() - 1;
					gWireItr2 = gnet._gWires.end() - 1;
				}
				gnet._gWires.push_back(
						GlobalWire(g[vecPath.back()].x, g[vecPath.back()].y,
								g[vecPath.back()].z,
								gWireItr1 - gnet._gWires.begin(), -1));
				gWireItr2->_pWireId = gnet._gWires.size() - 1;
			}
		} else {
// cout << "Mode 1" << endl;
// Mode 1
			if (vecPath.size() == 1) {
				gWireItr1->_pWireId = gWireItr2 - gnet._gWires.begin();
			} else {
				std::reverse(vecPath.begin(), vecPath.end());
				for (auto itr = vecPath.begin(); itr < vecPath.end() - 2;
						++itr) {
// Add edge (*itr, *(itr+1)) to DB (gnet/addEdgeDemand)
					gnet._gWires.push_back(
							GlobalWire(g[*itr].x, g[*itr].y, g[*itr].z, -1,
									-1));
					gWireItr1->_pWireId = gnet._gWires.size() - 1;
					gWireItr1 = gnet._gWires.end() - 1;
				}
				gnet._gWires.push_back(
						GlobalWire(g[vecPath[vecPath.size() - 2]].x,
								g[vecPath[vecPath.size() - 2]].y,
								g[vecPath[vecPath.size() - 2]].z,
								gWireItr2 - gnet._gWires.begin(), -1));
				gWireItr1->_pWireId = gnet._gWires.size() - 1;
			}
		}
		restoreNet(vp1.gnetID, gnet, _gnets[vp1.gnetID]._cleanEdges,
				_gnets[vp1.gnetID]._cleanVias, false);
	}
	return make_tuple(wl + 1, addedEdges, addedVias, margin); // including WL of the new v-pin
}

