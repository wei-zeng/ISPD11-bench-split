/*
 * outputCSV.cpp
 *
 *  Created on: Aug 25, 2019
 *      Author: wzeng
 */

#include "grDB.h"
#include "graph.h" // for undirected Graph from boost
#include <iterator>
#include <fstream>
#include <cassert>
#include <set>
#include <boost/graph/connected_components.hpp>

vector<Vpin> RoutingDB::getVpins(const Layout &layout, int splitLayer,
		bool twoCutNetsOnly, bool includeFloatingVpins, bool excludeNI) const { // splitLayer: 0-based
	vector<Vpin> vpins;
	using namespace boost;
	int totalVpin = 0;
	cout << "Splitting... " << endl;
	for (size_t gnetID = 0; gnetID < _gnets.size(); ++gnetID) {
		auto &gnet = _gnets[gnetID];
		Graph routingTree;
		for (auto &gWire : gnet._gWires) {
			Vertex v = add_vertex(routingTree);
			routingTree[v].x = gWire._x;
			routingTree[v].y = gWire._y;
			routingTree[v].z = gWire._z;
			routingTree[v].realPinId = gWire._realPinId;
			routingTree[v].isPin = (!excludeNI || gWire._z == 0) && (gWire._realPinId >= 0);
			routingTree[v].isRoot = (gWire._pWireId == -1
					&& gWire._realPinId >= 0);
		}
		for (size_t i = 0; i < gnet._gWires.size(); i++) {
			auto &gWire = gnet._gWires[i];
			if (gWire._pWireId >= 0
					&& !(gWire._z > splitLayer
							&& gnet._gWires[gWire._pWireId]._z > splitLayer)) {
				add_edge(vertex(gWire._pWireId, routingTree),
						vertex(i, routingTree), routingTree); // only connect if lower than or on splitLayer
			}
		}
		int vpCount = 0;
		EdgeIter e, eend;
		for (tie(e, eend) = edges(routingTree); e != eend; ++e) { // for each vpin
			if (routingTree[e->m_source].z > splitLayer
					|| routingTree[e->m_target].z > splitLayer) {
				routingTree[*e].isVpin = true;
				++vpCount;
				++totalVpin;
			}
		}
		if (twoCutNetsOnly && vpCount != 2)
			continue;
		int vpCountBase = totalVpin - vpCount;

		vpCount = 0;
		for (tie(e, eend) = edges(routingTree); e != eend; ++e) { // for each vpin
			if (routingTree[*e].isVpin) {
				Vpin vp;
				vp.vPinID = vpCountBase + vpCount;
				vp.gnetID = gnetID;
				vp.xCoord = routingTree[e->m_source].x;
				vp.yCoord = routingTree[e->m_source].y;
				vp.zCoord = splitLayer;
				std::vector<int> c(num_vertices(routingTree));
				connected_components(routingTree, &c[0]);
				VertexIter v, vend;
				EdgeIter e2, e2end;
				for (tie(e2, e2end) = edges(routingTree); e2 != e2end; ++e2) {
					if (c[e2->m_source] == c[e->m_source]) {
						int x1 = routingTree[e2->m_source].x;
						int y1 = routingTree[e2->m_source].y;
						int z1 = routingTree[e2->m_source].z;
						int x2 = routingTree[e2->m_target].x;
						int y2 = routingTree[e2->m_target].y;
						int z2 = routingTree[e2->m_target].z;
						if (x1 != x2) {
							assert(y1 == y2);
							assert(z1 == z2);
							for (int x = min(x1, x2); x <= max(x1, x2); ++x) {
								vp.gnetGrids.insert(Gcell(x, y1, z1));
							}
                            vp.routes.push_back(make_tuple(x1, y1, z1, x2, y2, z2));
						} else if (y1 != y2) {
							assert(x1 == x2);
							assert(z1 == z2);
							for (int y = min(y1, y2); y <= max(y1, y2); ++y) {
								vp.gnetGrids.insert(Gcell(x1, y, z1));
							}
                            vp.routes.push_back(make_tuple(x1, y1, z1, x2, y2, z2));
						} else if (z1 != z2) {
							assert(x1 == x2);
							assert(y1 == y2);
							for (int z = min(z1, z2); z <= max(z1, z2); ++z) {
								vp.gnetGrids.insert(Gcell(x1, y1, z));
							}
                            vp.routes.push_back(make_tuple(x1, y1, z1, x2, y2, z2));
						}
					}
				}
				for (tie(v, vend) = vertices(routingTree); v != vend; ++v) {
					if (c[*v] == c[e->m_source]) {
						vp.gnetVertices.insert(
								Gcell(routingTree[*v].x, routingTree[*v].y,
										routingTree[*v].z));
						vp.gnetGrids.insert(
								Gcell(routingTree[*v].x, routingTree[*v].y,
										routingTree[*v].z));
						if (routingTree[*v].isRoot)
							vp.hasOriginalRoot = true;
					}

					if (c[*v] == c[e->m_source]
							&& routingTree[*v].isPin) {
						Gcell gcell = Gcell(routingTree[*v].x,
								routingTree[*v].y, routingTree[*v].z);
						// find pin type
						// find pin node
						Net net = gnet._net;
						assert(net._netPins.size() == net._netPinTypes.size());
						assert(
								net._netPins.size() == net._netPinNodeIDs.size());
						for (size_t i = 0; i < net._netPins.size(); ++i) {
							if (layout.getGcell(net._netPins[i]) == gcell) {
								vp.pinNodeIDs.push_back(net._netPinNodeIDs[i]);
								vp.pins.push_back(net._netPins[i]);
								vp.pinTypes.push_back(net._netPinTypes[i]);
								auto &the_node =
										layout._nodes[net._netPinNodeIDs[i]];
								if (the_node._isTerm
										&& the_node._numInputs == 1
										&& the_node._numOutputs == 0) { // PI
									assert(vp.pinTypes.back() == CI);
									vp.pinTypes.back() = PI;
								} /*else if (the_node._isTerm
										&& the_node._numInputs == 0
										&& the_node._numOutputs == 1) { // PO
									assert(vp.pinTypes.back() == CO);
									vp.pinTypes.back() = PO;
								}*/
								vp.cellAreas.push_back(
										layout._nodes[net._netPinNodeIDs[i]]._height
												* layout._nodes[net._netPinNodeIDs[i]]._width);
							}
						}
					}
				}
				EdgeIter ewl, ewlend;
				int wl = 0;
				for (tie(ewl, ewlend) = edges(routingTree); ewl != ewlend;
						++ewl) {
					if (c[ewl->m_source] == c[e->m_source]
							&& c[ewl->m_target] == c[e->m_source]) {
						wl += Gcell(routingTree[ewl->m_source].x,
								routingTree[ewl->m_source].y,
								routingTree[ewl->m_source].z)
								- Gcell(routingTree[ewl->m_target].x,
										routingTree[ewl->m_target].y,
										routingTree[ewl->m_target].z);
					}
				}
				if (includeFloatingVpins || vp.pins.size() > 0) {
					vpCount++;
					vp.wlToL1 = wl;
					vpins.push_back(vp);
				}
			}
		}
		if (twoCutNetsOnly && vpCount != 2) {
			while (vpCount--) {
				vpins.pop_back();
			}
		}
		// Find matching Vpin
		for (int i = 1; i <= vpCount; ++i) {
			for (int j = 1; j <= vpCount; ++j) {
				if (vpins[vpins.size() - i].vPinID
						!= vpins[vpins.size() - j].vPinID) {
					vpins[vpins.size() - i].matchingVpinIdx = vpins.size() - j;
					break;
				}
			}
		}
	}

	cout << "# Vpins: " << vpins.size() << endl;
	return vpins;
}

