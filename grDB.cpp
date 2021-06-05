#include "grDB.h"
#include <unordered_set>
#include <random>
#include <algorithm>
#include <iterator>
#include <boost/graph/connected_components.hpp>
#include <iomanip>

//  Read Global Routes and initialize _gnets data structure
void RoutingDB::initGlobalNets(Layout &layout) {
	const int size = layout._nets.size();
	_gnets.resize(size);
	set<Gcell> gPins;
	set<Gcell> gPins2D;
	for (int i = 0; i < size; i++) {

		gPins.clear();
		gPins2D.clear();

		Net &net = layout._nets[i];
		for (auto pin : net._netPins) {
			Gcell gpin = layout.getGcell(pin);
			gPins2D.insert(Gcell(gpin._x, gpin._y, 0));
			gPins.insert(gpin);
		}

		_gnets[i]._name = net.getName();
		_gnets[i]._net = net;

		if (gPins2D.size() < 2)
			continue;

		vector<Gcell> pins(gPins.begin(), gPins.end());
		_gnets[i]._isGlobal = true;
		_gnets[i]._gPins.swap(pins);

		vector<Gcell> pins2D(gPins2D.begin(), gPins2D.end());
		_gnets[i]._gPins2D.swap(pins2D);

		++_numGlobalNets;
	}

	printf("<I> %-20s : %u\n", "# Global Nets", _numGlobalNets);
}

// Read routing file and initialize list of gnet points
void RoutingDB::readGlobalWires(const Layout &layout, const char *wiresFile) {
	readGlobalWires(layout, string(wiresFile));
}

// Read routing file and initialize list of gnet points
void RoutingDB::readGlobalWires(const Layout &layout, const string &wiresFile) {
	ifstream inp(wiresFile);
	if (!inp.is_open()) {
		printf("ERROR %s %s %d\n", __FILE__, __func__, __LINE__);
		exit(1);
	}

	initVias(layout);
	string line;
	int netId = 0;
	float x1, y1, x2, y2;
	int z1, z2;
	while (getline(inp, line)) {
		if (line == "!")
			continue;
		if (line.at(0) == 'n') {
			sscanf(line.c_str(), "n%*u %u", &netId);
			continue;
		}
		sscanf(line.c_str(), "(%f,%f,%d)-(%f,%f,%d)", &x1, &y1, &z1, &x2, &y2,
				&z2);

		if (z1 == z2) {
			if (x1 != x2)
				assert(y1 == y2);
			if (y1 != y2)
				assert(x1 == x2);
		} else {
			assert(x1 == x2 && y1 == y2);
		}
		_gnets[netId].addPtPair(layout.getGcell(Point3D(x1, y1, z1 - 1)),
				layout.getGcell(Point3D(x2, y2, z2 - 1)));
		_gnets[netId]._isRealNet = true;
	}
	inp.close();
}

void RoutingDB::writeNets(const Layout &layout, const char *netsFile, const vector<Vpin> &vpins, bool brokenNetsOnly) const {
	unordered_set<size_t> excludedGnetIds;
	for (auto vp : vpins) {
		excludedGnetIds.insert(vp.gnetID);
	}
	ofstream ofs(netsFile);
	cout << "Writing nets" << endl;
	size_t netCount = 0, pinCount = 0;
	if (!brokenNetsOnly) {
		for (size_t i = 0; i < layout._nets.size(); ++i) {
			if (excludedGnetIds.count(i))
				continue;
			netCount++;
			pinCount += layout._nets.at(i)._netPins.size();
		}
	}
	for (auto &vp : vpins) {
		netCount++;
		pinCount += vp.pins.size();
	}
	ofs << "UCLA nets 1.0\n\n" << "NumNets  :  " << netCount << endl
			<< "NumPins  :  " << pinCount << endl << endl;
	if (!brokenNetsOnly) {
		for (size_t i = 0; i < layout._nets.size(); ++i) {
			if (excludedGnetIds.count(i))
				continue;
			for (auto line : layout._nets.at(i).infoLines()) {
				ofs << line << endl;
			}
		}
	}
	for (auto &vp : vpins) {
		ofs << "NetDegree  :  " << vp.pins.size() << "    n" << vp.gnetID << "_" << vp.vPinID << endl;
		for (size_t i = 0; i < vp.pinNodeIDs.size(); ++i) {
			const Node &node = layout._nodes.at(vp.pinNodeIDs.at(i));
			ofs << setw(15) << node._name << "   " << (vp.pinTypes[i] == CI || vp.pinTypes[i] == PI ? "I" : vp.pinTypes[i] == CO || vp.pinTypes[i] == PO ? "O" : "N") << "  :";
			ofs << " " << setw(11) << fixed << setprecision(4) << 
				vp.pins[i]._x - node.getLLX() - node.getWidth() / 2;
			ofs << " " << setw(11) << fixed << setprecision(4) << 
				vp.pins[i]._y - node.getLLY() - node.getHeight() / 2 << endl;
		}
	}
}

void RoutingDB::writeKey(const Layout &layout, const char *keyFile, const vector<Vpin> &vpins) const {
	unordered_set<size_t> excludedGnetIds;
	for (auto vp : vpins) {
		excludedGnetIds.insert(vp.gnetID);
	}
	ofstream ofs(keyFile);
	cout << "Writing key" << endl;
	ofs << "SplitMan Key\n" << endl;
	for (auto &vp : vpins) {
		ofs << "n" << vp.gnetID << "_" << vp.vPinID << "  :  "
			<< "n" << vp.gnetID << endl;
	}
}


size_t RoutingDB::writeGlobalWires(const Layout &layout,
		const char *wiresFile) const {
    unordered_set<size_t> empty_set;
	return writeGlobalWires(layout, wiresFile, empty_set);
}

// Read routing file and initialize list of gnet points
size_t RoutingDB::writeGlobalWires(const Layout &layout,
		const char* wiresFile, const unordered_set<size_t> &excludedGnetIds) const {
	string wiresFile_str = string(wiresFile);
    ofstream outp(wiresFile_str);
	if (!outp.is_open()) {
		printf("ERROR %s %s %d\n", __FILE__, __func__, __LINE__);
		exit(1);
	}
	int xGrid = layout._cellWidth, yGrid = layout._cellHeight;
	char buf[1024];
	float x1, y1, x2, y2;
	int z1, z2;
    size_t netCount = 0;
	for (size_t netId = 0; netId < _gnets.size(); netId++) {
		if (!_gnets[netId]._isRealNet)
			continue;
        if (excludedGnetIds.count(netId))
            continue;
		auto &gnet = _gnets[netId];
			sprintf(buf, "n%lu %lu", netId, netCount++);
			outp << buf << endl;
			for (auto gw : gnet._gWires) {
				if (gw._pWireId != -1) {
					x1 = gw._x;
					y1 = gw._y;
					z1 = gw._z;
					x2 = gnet._gWires[gw._pWireId]._x;
					y2 = gnet._gWires[gw._pWireId]._y;
					z2 = gnet._gWires[gw._pWireId]._z;
					if (x1 > x2 || y1 > y2 || z1 > z2) {
						swap(x1, x2);
						swap(y1, y2);
						swap(z1, z2);
					}
					sprintf(buf, "(%.0f,%.0f,%d)-(%.0f,%.0f,%d)", x1 * xGrid,
							y1 * yGrid, z1 + 1, x2 * xGrid, y2 * yGrid, z2 + 1);
					outp << buf << endl;
				}
			}
			//sprintf(buf, "(%.0f,%.0f,%d)-(%.0f,%.0f,%d)", xs * xGrid,
			//		ys * yGrid, zs + 1, xf * xGrid, yf * yGrid, zf + 1);
			//outp << buf << endl;
			outp << "!" << endl;
		}
	outp.close();
    return netCount;
}

void RoutingDB::writeGlobalWiresVpin(const Layout &layout,
		const char* wiresFile, const Vpin &vp, size_t &netCount) const {
	string wiresFile_str = string(wiresFile);
    ofstream outp(wiresFile_str, ios::out | ios::app);
	if (!outp.is_open()) {
		printf("ERROR %s %s %d\n", __FILE__, __func__, __LINE__);
		exit(1);
	}
	int xGrid = layout._cellWidth, yGrid = layout._cellHeight;
	char buf[1024];
	int x1_int, y1_int, x2_int, y2_int;
    float x1, y1, x2, y2;
	int z1, z2;
	sprintf(buf, "n%d_%d %lu", vp.gnetID, vp.vPinID, netCount++);
	outp << buf << endl;
	for (auto rt : vp.routes) {
		tie(x1_int, y1_int, z1, x2_int, y2_int, z2) = rt;
		x1 = x1_int;
		y1 = y1_int;
		x2 = x2_int;
		y2 = y2_int;
		if (x1 > x2 || y1 > y2 || z1 > z2) {
			swap(x1, x2);
			swap(y1, y2);
			swap(z1, z2);
		}
		sprintf(buf, "(%.0f,%.0f,%d)-(%.0f,%.0f,%d)", x1 * xGrid,
				y1 * yGrid, z1 + 1, x2 * xGrid, y2 * yGrid, z2 + 1);
		outp << buf << endl;
	}
	outp << "!" << endl;
	outp.close();
}


// Initialize routing tree for every Gnet
void RoutingDB::initRoutingTreeForAllGnets() {
	_treeTable.resize(_numLayers);
	_uMapCell2LocInTable.resize(_numLayers);

	for (size_t i = 0; i < _gnets.size(); i++) {
		Gnet &gnet = _gnets[i];
		if (!gnet.isGlobal())
			continue;
		initRoutingTree(gnet);
	}
	_treeTable.clear();
	_uMapCell2LocInTable.clear();
}

// Initialize routing tree, populate treeTable for routing purposes
// Structure of treeTable : (int x, int y, int z, int indexId, int distPrev, int distNext, bool hasDnVia, bool hasUpVia)
void RoutingDB::initRoutingTree(Gnet &gnet) {
	for (int z = 0; z < _numLayers; z++) {
		_treeTable[z].clear();
		_uMapCell2LocInTable[z].clear();
	}

	int indexId = 0;
	for (auto gPin : gnet._gPins) {
		_treeTable[gPin._z].push_back(
				PtInfo(gPin._x, gPin._y, gPin._z, indexId++, -1, -1, false,
						false));
	}

// Add points to treeTable
	for (auto ptPair : gnet._gpPairArray) {
		Gcell &pt1 = ptPair.first;
		Gcell &pt2 = ptPair.second;
		if (pt1 == pt2)
			continue;
		const int dist = pt1 - pt2;
		if (pt1._z == pt2._z) {
			_treeTable[pt1._z].push_back(
					PtInfo(pt1._x, pt1._y, pt1._z, -1, pt1 < pt2 ? -1 : dist,
							pt1 < pt2 ? dist : -1, false, false));
			_treeTable[pt2._z].push_back(
					PtInfo(pt2._x, pt2._y, pt2._z, -1, pt1 < pt2 ? dist : -1,
							pt1 < pt2 ? -1 : dist, false, false));
		} else {
			const int zMin = min(pt1._z, pt2._z);
			const int zMax = max(pt1._z, pt2._z);
			for (int z = zMin; z <= zMax; z++)
				_treeTable[z].push_back(
						PtInfo(pt1._x, pt1._y, z, -1, -1, -1, z != zMin,
								z != zMax));
		}
	}

	for (int z = 0; z < _numLayers; z++) {
		if (_treeTable[z].empty())
			continue;
		tmp.clear();
		if (_dirLayers[z] == V)
			sort(_treeTable[z].begin(), _treeTable[z].end(),
					comparePtInfoOnVerLayer());
		else
			sort(_treeTable[z].begin(), _treeTable[z].end(),
					comparePtInfoOnHorLayer());

		int ptNum = 0;
		for (vector<PtInfo>::iterator itr = _treeTable[z].begin();
				itr != _treeTable[z].end(); itr++) {
			if (tmp.empty() || tmp.back() != (*itr)) {
				_uMapCell2LocInTable[z].insert(
						make_pair(itr->getLoc(), ptNum++));
				tmp.push_back(
						PtInfo(itr->_x, itr->_y, itr->_z, -1, -1, -1, false,
								false));
			}
			PtInfo &pt = tmp.back();
			if (itr->_indexId != -1) {
				try {
					if (pt._indexId != -1)
						throw runtime_error("gPins cannot overlap in gNet");
				} catch (const runtime_error &e) {
					cout << "Exception Caught : " << e.what() << endl;
				}
				pt._indexId = itr->_indexId;
			}
			if (itr->_distPrev > pt._distPrev)
				pt._distPrev = itr->_distPrev;
			if (itr->_distNext > pt._distNext)
				pt._distNext = itr->_distNext;
			if (itr->_hasDnVia)
				pt._hasDnVia = true;
			if (itr->_hasUpVia)
				pt._hasUpVia = true;
		}

		_treeTable[z].swap(tmp);

		if (_treeTable[z].size() == 1)
			continue;

		for (int cur = 0; cur < static_cast<int>(_treeTable[z].size()); cur++) {
			PtInfo &pt = _treeTable[z][cur];
			if (pt._distNext > -1) {
				int next = _uMapCell2LocInTable[z][
						_dirLayers[z] == H ?
								Gcell(pt._x + pt._distNext, pt._y, pt._z) :
								Gcell(pt._x, pt._y + pt._distNext, pt._z)];
				for (int i = cur; i < next; i++) {
					int dist = _treeTable[z][i] - _treeTable[z][i + 1];
					_treeTable[z][i]._distNext = dist;
					_treeTable[z][i + 1]._distPrev = dist;
				}
			}
		}
	}

//  add root gwire
	gnet._gWires.push_back(
			GlobalWire(gnet._gPins[0]._x, gnet._gPins[0]._y, gnet._gPins[0]._z,
					-1, 0));
	gnet._root = 0;
	frontier.clear();
	frontier.push_back(gnet._gPins[0]);
	size_t head = 0;
	while (head < frontier.size()) {
		Gcell &pt = frontier[head];
		int loc = _uMapCell2LocInTable[pt._z][pt];
		PtInfo &ptInfo = _treeTable[pt._z][loc];
		if (ptInfo._distNext > -1) {
			PtInfo &neighorPtInfo = _treeTable[pt._z][loc + 1];
			neighorPtInfo._distPrev = -1;
			gnet._gWires.push_back(
					GlobalWire(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z, head, neighorPtInfo._indexId));
			frontier.push_back(
					Gcell(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z));
		}
		if (ptInfo._distPrev > -1) {
			PtInfo &neighorPtInfo = _treeTable[pt._z][loc - 1];
			neighorPtInfo._distNext = -1;
			gnet._gWires.push_back(
					GlobalWire(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z, head, neighorPtInfo._indexId));
			frontier.push_back(
					Gcell(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z));
		}
		if (ptInfo._hasUpVia) {
			int otherLoc = _uMapCell2LocInTable[pt._z + 1][Gcell(pt._x, pt._y,
					pt._z + 1)];
			PtInfo &neighorPtInfo = _treeTable[pt._z + 1][otherLoc];
			neighorPtInfo._hasDnVia = false;
			gnet._gWires.push_back(
					GlobalWire(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z, head, neighorPtInfo._indexId));
			frontier.push_back(
					Gcell(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z));
		}
		if (ptInfo._hasDnVia) {
			int otherLoc = _uMapCell2LocInTable[pt._z - 1][Gcell(pt._x, pt._y,
					pt._z - 1)];
			PtInfo &neighorPtInfo = _treeTable[pt._z - 1][otherLoc];
			neighorPtInfo._hasUpVia = false;
			gnet._gWires.push_back(
					GlobalWire(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z, head, neighorPtInfo._indexId));
			frontier.push_back(
					Gcell(neighorPtInfo._x, neighorPtInfo._y,
							neighorPtInfo._z));
		}
		head++;
	}
#if 0
    /* read through the wire vector and remove antenna wires in linear time */
    _numChildren.assign( gnet._gWires.size(), 0 );
    _offset.assign( gnet._gWires.size(), 0 );
    _validWires.assign( gnet._gWires.size(), true );
    for( size_t i = 1; i < gnet._gWires.size(); i++ )
        _numChildren[gnet._gWires[i]._pWireId] ++;
    for( size_t i = 0; i < gnet._gWires.size(); i++ ) {
        if( _numChildren[i] == 0 && gnet._gWires[i]._realPinId == -1 ) {
            /* antenna wire found */
            int curr = i;
            while( curr != -1 && _numChildren[curr] == 0 && gnet._gWires[curr]._realPinId == -1 && _validWires[curr] ) {
                _validWires[curr] = false;
                curr = gnet._gWires[curr]._pWireId;
                if( curr != -1 )
                    _numChildren[curr] --;
            }
        }
    }
    int currOffset = 0;
    int curr = 0;
    _cleanWires.resize( gnet._gWires.size() );
    for( size_t i = 0; i < gnet._gWires.size(); i++ ) {
    	if( _validWires[i] ) {
            _cleanWires[curr] = gnet._gWires[i];
            if( _cleanWires[curr]._pWireId != -1 ) {
                _cleanWires[curr]._pWireId -= _offset[ _cleanWires[curr]._pWireId ];
            }
            curr ++;
        }
        else {
            currOffset ++;
        }
        _offset[i] = currOffset;
    }
    assert( curr + currOffset == static_cast<int>(gnet._gWires.size()) );
    _cleanWires.resize( curr );
    gnet._gWires.swap( _cleanWires );
#endif
}

// Update edge and via demands from info on _gWires
void RoutingDB::updateEdgeDemands() {
	for (size_t k = 0; k < _gnets.size(); ++k) {
		auto &gnet = _gnets[k];
		if (!gnet.isGlobal())
			continue;
		for (size_t i = 1; i < gnet._gWires.size(); i++) {
			const GlobalWire &gWire = gnet._gWires[i];
			const GlobalWire &pWire = gnet._gWires[gWire._pWireId];
			if (pWire._z == gWire._z) {
				if (_dirLayers[pWire._z] == H) {
					assert(gWire._y == pWire._y);
					for (int x = min(gWire._x, pWire._x);
							x < max(gWire._x, pWire._x); x++) {
						size_t edgeId = findEdge(x, gWire._y, gWire._z);
						addEdgeDemand(edgeId);
						gnet._occupiedEdges.insert(edgeId);
					}
				} else if (_dirLayers[pWire._z] == V) {
					assert(gWire._x == pWire._x);
					for (int y = min(gWire._y, pWire._y);
							y < max(gWire._y, pWire._y); y++) {
						size_t edgeId = findEdge(gWire._x, y, gWire._z);
						addEdgeDemand(edgeId);
						gnet._occupiedEdges.insert(edgeId);
					}
				} else
					assert(0);
			} else {
				assert(gWire._y == pWire._y);
				assert(gWire._x == pWire._x);
				for (int z = min(gWire._z, pWire._z);
						z < max(gWire._z, pWire._z); z++) {
					size_t viaId = findVia(gWire._x, gWire._y, z);
					addViaDemand(viaId);
					gnet._occupiedVias.insert(viaId);
				}
			}
		}
	}
}

// Create routing blockage and store info in layout data structure
void RoutingDB::createBlockage(const int xStart, const int yStart,
		const int xEnd, const int yEnd, const int z, Layout &layout) {

	int availableEdgeCap = 0;

	if (_dirLayers[z] == H) {
		for (int y = yStart; y < yEnd; y++) {
			for (int x = xStart; x < xEnd; x++) {
				const int eid = findEdge(x, y, z);
				availableEdgeCap += max(0,
						getEdgeCap(eid) - getEdgeBlk(eid) - getEdgeDemand(eid));
			}
		}
	} else if (_dirLayers[z] == V) {
		for (int x = xStart; x < xEnd; x++) {
			for (int y = yStart; y < yEnd; y++) {
				const int eid = findEdge(x, y, z);
				availableEdgeCap += max(0,
						getEdgeCap(eid) - getEdgeBlk(eid) - getEdgeDemand(eid));
			}
		}
	}

	int size = 0;
	while (size * size * _capOnLayer[z] < availableEdgeCap) {
		size++;
	}

	char objName[256];
	sprintf(objName, "o%d", layout._numRegularNodes);
	layout._numRegularNodes++;

	layout.addNodeBlock(string(objName), size * _tileWidth, size * _tileHeight,
			xStart * _tileWidth, yStart * _tileHeight, z);
}

int RoutingDB::getTotalWL() const {
	return getWlByNet(-1);
}

int RoutingDB::getWlByNet(int gnetID) const { // gnetID < 0 for all gnets
	int wl = 0;
	int via = 0;
	size_t start = (gnetID < 0 ? 0 : gnetID);
	size_t end = (gnetID < 0 ? _gnets.size() : gnetID + 1);
	for (size_t i = start; i < end; i++) {
		auto &gnet = _gnets[i];
		if (!gnet.isGlobal())
			continue;
		for (size_t i = 1; i < gnet._gWires.size(); i++) {
			const GlobalWire &gWire = gnet._gWires[i];
			if (gWire._pWireId == -1)
				continue;
			const GlobalWire &pWire = gnet._gWires[gWire._pWireId];
			if (pWire._z == gWire._z) {
				if (_dirLayers[pWire._z] == H) {
					assert(gWire._y == pWire._y);
					for (int x = min(gWire._x, pWire._x);
							x < max(gWire._x, pWire._x); x++) {
						wl++;
					}
				} else if (_dirLayers[pWire._z] == V) {
					assert(gWire._x == pWire._x);
					for (int y = min(gWire._y, pWire._y);
							y < max(gWire._y, pWire._y); y++) {
						wl++;
					}
				} else
					assert(0);
			} else {
				assert(gWire._y == pWire._y);
				assert(gWire._x == pWire._x);
				for (int z = min(gWire._z, pWire._z);
						z < max(gWire._z, pWire._z); z++) {
					via++;
				}
			}
		}
	}
	printf("<I> Wirelength: %d, Via: %d, Total: %d\n", wl, via, wl + via);
	cout << flush;
	return wl + via;
}

int RoutingDB::getTotalWLByDemand() const {
	int wl = 0;
	int via = 0;
	for (size_t e = 0; e < _edges.size(); e++) {
		wl += _edges[e].demand() / _trackDemands[getEdgeLayer(e)];
	}
	for (size_t v = 0; v < _vias.size(); v++) {
		via += _vias[v].demand() / _trackDemands[getViaLayer(v) + 1]
				/ _trackDemands[getViaLayer(v) + 1];
	}
	printf("<I> From demand: Wirelength: %d, Via: %d, Total: %d\n", wl, via,
			wl + via);
	cout << flush;
	return wl + via;
}
