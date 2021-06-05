#include "grDB.h"

bool congestionEqual(const RoutingDB &db1, const RoutingDB &db2) {
	auto es1 = db1.getEdges();
	auto es2 = db2.getEdges();
	if (es1.size() != es2.size()) {
		cerr << "Via size different!" << endl;
		return false;
	}
	for (size_t i = 0; i < es1.size(); i++) {
		if (!(es1[i] == es2[i])) {
			cerr << i << "-th edge different" << endl;
			Gcell gc = db1.getEdgeLoc(i);
			cerr << gc._x << "," << gc._y << "," << gc._z << endl;
			cerr << "Demand: old " << es1[i].demand() << " vs new "
					<< es2[i].demand() << endl;
			cerr << "Cap: old " << es1[i].cap() << " vs new " << es2[i].cap()
					<< endl;
			cerr << "Block: old " << es1[i].blk() << " vs new " << es2[i].blk()
					<< endl;
			return false;
		}
	}
	auto vs1 = db1.getVias();
	auto vs2 = db2.getVias();
	if (vs1.size() != vs2.size()) {
		cerr << "Via size different!" << endl;
		return false;
	}
	for (size_t i = 0; i < vs1.size(); i++) {
		if (!(vs1[i] == vs2[i])) {
			cerr << i << "-th via different" << endl;
			return false;
		}
	}
	return true; //db1.getEdges() == db2.getEdges() && db1.getVias() == db2.getVias();
}

//  
//  Edge Related Functions
//
void RoutingDB::initEdges(const Layout &layout) {
	_edges.resize(
			max(_numTilesX * (_numTilesY - 1) * _numLayers,
					(_numTilesX - 1) * _numTilesY * _numLayers));
	_cntEdgeLayers.assign(_numLayers + 1, 0);
	_trackDemands.assign(_numLayers, 0);
	_dirLayers.assign(_numLayers, NA);
	_capOnLayer.assign(_numLayers, 0);
	_upViaSize.assign(_numLayers + 1, 0);

	//  non-via
	int cnt = 0;
	for (int layer = 0; layer < _numLayers; layer++) {
		// routing dir on that layer
		_dirLayers[layer] = layout._hCaps[layer] > 0 ? H :
							layout._vCaps[layer] > 0 ? V : NA;
		//  track demand on that layer
		_trackDemands[layer] = static_cast<int>(layout._mww[layer]
				+ layout._mws[layer]);
		//  edge count marks
		_cntEdgeLayers[layer + 1] = _cntEdgeLayers[layer]
				+ (_dirLayers[layer] == H ? _numTilesY * (_numTilesX - 1) :
					_dirLayers[layer] == V ? _numTilesX * (_numTilesY - 1) : 0);

		//  via size
		if (layer > 0) {
			_upViaSize[layer - 1] = _trackDemands[layer];
		}

		//  initialize the free lists
		if (_dirLayers[layer] == H) {
			for (int y = 0; y < _numTilesY; y++)
				for (int x = 0; x < _numTilesX - 1; x++)
					_edges[cnt++].initFreeList(y * layout._cellHeight,
							(y + 1) * layout._cellHeight);
		}
		if (_dirLayers[layer] == V) {
			for (int x = 0; x < _numTilesX; x++)
				for (int y = 0; y < _numTilesY - 1; y++)
					_edges[cnt++].initFreeList(x * layout._cellWidth,
							(x + 1) * layout._cellWidth);
		}
	}
	_upViaSize[_numLayers] = _upViaSize[_numLayers - 1];
	printf("<I> %-20s : ", "Demand Per Track");
	for (int i = 0; i < _numLayers; i++)
		printf("%6d", _trackDemands[i]);
	printf("\n");
	_edges.resize(_cntEdgeLayers[_numLayers]);
}

RoutingDir RoutingDB::getRoutingDir(const int layer) const {
	return _dirLayers[layer];
}

int RoutingDB::getEdgeLayer(const int edgeId) const {
	int layer = 0;
	while (_cntEdgeLayers[layer + 1] <= edgeId)
		layer++;
	return layer;
}

int RoutingDB::getEdgeDemand(const int edgeId) const {
	return _edges[edgeId].demand();
}

int RoutingDB::getEdgeBlk(const int edgeId) const {
	return _edges[edgeId].blk();
}

int RoutingDB::getEdgeCap(const int edgeId) const {
	return _edges[edgeId].cap();
}

int RoutingDB::getEdgeOf(const int edgeId) const {
	return _edges[edgeId].of();
}

Gcell RoutingDB::getEdgeLoc(const int edgeId) const {
	int layer = 0;
	while (_cntEdgeLayers[layer + 1] <= edgeId)
		layer++;
	int edgeIdOnLayer = edgeId - _cntEdgeLayers[layer];
	if (_dirLayers[layer] == H) {
		return Gcell(edgeIdOnLayer % (_numTilesX - 1),
				edgeIdOnLayer / (_numTilesX - 1), layer);
	} else {
		return Gcell(edgeIdOnLayer / (_numTilesY - 1),
				edgeIdOnLayer % (_numTilesY - 1), layer);
	}
}

void RoutingDB::addEdgeDemand(const int edgeId) {
	_edges[edgeId].addDemand(_trackDemands[getEdgeLayer(edgeId)]);
}

void RoutingDB::reduceEdgeDemand(const int edgeId) {
	_edges[edgeId].addDemand(-_trackDemands[getEdgeLayer(edgeId)]);
}

void RoutingDB::ripUpEdge(const int edgeId) {
	reduceEdgeDemand(edgeId);
}

void RoutingDB::restoreEdge(const int edgeId) {
	addEdgeDemand(edgeId);
}

// Set _blk for all edges
void RoutingDB::initGlobalEdgeProfile(const Layout &layout) {
	for (size_t i = 0; i < layout._blkId.size(); i++) {
		const Node &blkNode = layout._nodes[layout._blkId[i]];
		refreshNodeBlkInfo(layout, blkNode);
	}
	int edgeId = 0;
	for (int layer = 0; layer < _numLayers; layer++) {
		_capOnLayer[layer] =
				_dirLayers[layer] == H ? layout._hCaps[layer] :
				_dirLayers[layer] == V ? layout._vCaps[layer] : 0;
		if (_dirLayers[layer] == NA)
			continue;
		else {
			const double fullLength =
					_dirLayers[layer] == H ?
							layout._cellHeight : layout._cellWidth;
			const int edgeNumOnLayer =
					_dirLayers[layer] == H ?
							_numTilesY * (_numTilesX - 1) :
							_numTilesX * (_numTilesY - 1);
			for (int i = 0; i < edgeNumOnLayer; i++) {
				Edge &edge = _edges[edgeId];
				edge.setCap(_capOnLayer[layer]);
				int numAvailableTracks = floor(
						edge.cap() * (edge.getFreeLength() / fullLength)
								/ _trackDemands[layer]);
				edge.setBlk(
						edge.cap() - _trackDemands[layer] * numAvailableTracks);
				edgeId++;
			}
		}
	}
	printf("<I> %-20s : %d\n", "# Global Edges", edgeId);
}

int RoutingDB::findEdge(const int x, const int y, const int z) const {
	return _dirLayers[z] == H ?
			_cntEdgeLayers[z] + y * (_numTilesX - 1) + x :
			_cntEdgeLayers[z] + x * (_numTilesY - 1) + y;
}

int RoutingDB::getViaLayer(const int viaId) const {
	return viaId / (_numTilesX * _numTilesY);
}

int RoutingDB::getViaDemand(const int viaId) const {
	return _vias[viaId].demand();
}

int RoutingDB::getViaCap(const int viaId) const {
	return _vias[viaId].cap();
}

int RoutingDB::getViaOf(const int viaId) const {
	return _vias[viaId].of();
}

Gcell RoutingDB::getViaLoc(const int viaId) const {
	return Gcell(viaId % _numTilesX, (viaId / _numTilesX) % _numTilesY,
			viaId / (_numTilesX * _numTilesY));
}

//
//  Blockage Related Functions
//
void RoutingDB::refreshNodeBlkInfo(const Layout &layout, const Node &blkNode) {
	if (!blkNode.isNonRect()) {
		const Shape rect(blkNode.getLLX(), blkNode.getLLY(), blkNode.getWidth(),
				blkNode.getHeight());
		for (size_t i = 0; i < blkNode._blockedLayers.size(); i++)
			refreshRectBlkInfo(layout, rect, blkNode._blockedLayers[i]);
	} else {
		for (auto shape : blkNode._shapes)
			for (size_t i = 0; i < blkNode._blockedLayers.size(); i++)
				refreshRectBlkInfo(layout, shape, blkNode._blockedLayers[i]);
	}
}

void RoutingDB::refreshRectBlkInfo(const Layout &layout, const Shape &shape,
		const int layerId) {
	if (_dirLayers[layerId] == NA)
		return;

	const int xMin = static_cast<int>(shape._llx) / layout._cellWidth;
	const int yMin = static_cast<int>(shape._lly) / layout._cellHeight;
	const int xMax = static_cast<int>(shape._llx + shape._width)
			/ layout._cellWidth;
	const int yMax = static_cast<int>(shape._lly + shape._height)
			/ layout._cellHeight;

	for (int x = xMin; x <= xMax; x++)
		for (int y = yMin; y <= yMax; y++)
			updateEdgeFreeList(layout, findEdge(x, y, layerId), shape);
}

void RoutingDB::updateEdgeFreeList(const Layout &layout, const int edgeId,
		const Shape &shape) {
	const Gcell edgeLoc = getEdgeLoc(edgeId);
	double blkStart = 0;
	double blkEnd = 0;
	if (_dirLayers[edgeLoc._z] == H) {
		if ((edgeLoc._x + 1) * layout._cellWidth < shape._llx
				|| (edgeLoc._x + 1) * layout._cellWidth
						> shape._llx + shape._width || shape._height == 0)
			return;
		blkStart = shape._lly;
		blkEnd = shape._lly + shape._height;
	} else if (_dirLayers[edgeLoc._z] == V) {
		if ((edgeLoc._y + 1) * layout._cellHeight < shape._lly
				|| (edgeLoc._y + 1) * layout._cellHeight
						> shape._lly + shape._height || shape._width == 0)
			return;
		blkStart = shape._llx;
		blkEnd = shape._llx + shape._width;
	} else
		assert(0);
	updateSingleEdgeFreeList(_edges[edgeId], blkStart, blkEnd);
}

void RoutingDB::updateSingleEdgeFreeList(Edge &edge, const double blkStart,
		const double blkEnd) {
	list<pair<double, double> > &_freeList = edge._freeList;
	for (auto itr = _freeList.begin(); itr != _freeList.end();) {
		if (itr->second <= blkStart || itr->first >= blkEnd) {
			++itr;
		} else if (itr->first >= blkStart && itr->second <= blkEnd) {
			itr = _freeList.erase(itr);
		} else if (itr->first < blkStart && itr->second <= blkEnd) {
			itr->second = blkStart;
			++itr;
		} else if (itr->first >= blkStart && itr->second > blkEnd) {
			itr->first = blkEnd;
			++itr;
		} else {
			_freeList.insert(itr, make_pair(itr->first, blkStart));
			*itr = make_pair(blkEnd, itr->second);
			itr++;
		}
	}
}

void RoutingDB::checkEdgeFreeList(const int edgeId) {
	const Edge &edge = _edges[edgeId];
	const Gcell pt = getEdgeLoc(edgeId);
	const double start = static_cast<double>(
			_dirLayers[pt._z] == H ? pt._y * _tileHeight : pt._x * _tileWidth);
	const double end = start
			+ static_cast<double>(
					_dirLayers[pt._z] == V ? _tileWidth : _tileHeight);
	for (auto itr = edge._freeList.begin(); itr != edge._freeList.end();
			itr++) {
		try {
			if (itr->first > itr->second)
				throw runtime_error("invalid free list segment");
			if (itr->first < start)
				throw runtime_error("invalid free list segment head");
			if (itr->second > end)
				throw runtime_error("invalid free list segment tail");
		} catch (const runtime_error &e) {
			cout << "Exception Caught : " << e.what() << endl;
			printCompleteEdgeInfo(edgeId);
			exit(0);
		}
	}
}

// Prints all edge info including free list of edges
void RoutingDB::printCompleteEdgeInfo(const int edgeId) {
	const Edge &edge = _edges[edgeId];
	printf("<D> Loc (%4d,%4d,%4d) %c %8d%8d%8d%8.1f%20s\n",
			getEdgeLoc(edgeId)._x, getEdgeLoc(edgeId)._y, getEdgeLayer(edgeId),
			_dirLayers[getEdgeLayer(edgeId)] == H ? 'H' :
			_dirLayers[getEdgeLayer(edgeId)] == V ? 'V' : 'N', edge.cap(),
			edge.blk(), edge.demand(), edge.getFreeLength(),
			edge.blk() + edge.demand() > edge.cap() ? "OVERFLOW" : "");
	printf("<D> Print Edge Free List : \n");
	for (auto itr = edge._freeList.begin(); itr != edge._freeList.end();
			itr++) {
		printf("<D> \tFrom %f To %f\n", itr->first, itr->second);
	}
}

// Prints info on edges
void RoutingDB::printEdgeUsage(const size_t edgeId) {
	const Edge &edge = _edges[edgeId];
	printf("Edge( %d, %d, %d ) %c cap = %d, blk = %d, demand = %d, OF = %d",
			getEdgeLoc(edgeId)._x, getEdgeLoc(edgeId)._y, getEdgeLayer(edgeId),
			_dirLayers[getEdgeLayer(edgeId)] == H ? 'H' :
			_dirLayers[getEdgeLayer(edgeId)] == V ? 'V' : 'N', edge.cap(),
			edge.blk(), edge.demand(), edge.of());
}

// Prints info on edges
void RoutingDB::printViaUsage(const size_t viaId) {
	const Via &via = _vias[viaId];
	printf("Via( %d, %d, %d ) cap = %d, demand = %d, OF = %d",
			getViaLoc(viaId)._x, getViaLoc(viaId)._y, getViaLayer(viaId),
			via.cap(), via.demand(), via.of());
}

//  gWire related methods
bool RoutingDB::isUpViaGlobalWire(const Gnet &gnet,
		const GlobalWire &gWire) const {
	return (gnet._gWires[gWire._pWireId]._z > gWire._z);
}

bool RoutingDB::isDnViaGlobalWire(const Gnet &gnet,
		const GlobalWire &gWire) const {
	return (gnet._gWires[gWire._pWireId]._z < gWire._z);
}

void RoutingDB::bruteForceFindAllChildWires(const Gnet &gnet, const int wireId,
		vector<int> &cache) {
	cache.clear();
	for (int i = 0; i < static_cast<int>(gnet._gWires.size()); i++)
		if (gnet._gWires[i]._pWireId == wireId)
			cache.push_back(i);
}

void RoutingDB::clearEdgeDemands() {
	for (size_t i = 0; i < _edges.size(); i++)
		_edges[i].clearDemand();
}

//
//  Via Related Objects
//
void RoutingDB::initVias(const Layout &layout) {
	_vias.resize(_numTilesX * _numTilesY * (_numLayers - 1));
	int cnt = 0;
	for (int layer = 0; layer < _numLayers - 1; layer++) {
		const int cap =
				_dirLayers[layer + 1] == H ?
						layout._hCaps[layer + 1] : layout._vCaps[layer + 1];
		const int area = cap * cap;
		for (int j = 0; j < _numTilesX * _numTilesY; j++)
			_vias[cnt++].setCap(area);
	}
}

long RoutingDB::findVia(const int x, const int y, const int z) const {
	return z * _numTilesX * _numTilesY + y * _numTilesX + x;
}

void RoutingDB::addViaDemand(const int viaId) {
	int layer = 0;
	int cnt = _numTilesX * _numTilesY;
	while (cnt <= viaId) {
		cnt += _numTilesX * _numTilesY;
		layer++;
	}
	_vias[viaId].addDemand(_trackDemands[layer + 1] * _trackDemands[layer + 1]);
}

void RoutingDB::reduceViaDemand(const int viaId) {
	int layer = 0;
	int cnt = _numTilesX * _numTilesY;
	while (cnt <= viaId) {
		cnt += _numTilesX * _numTilesY;
		layer++;
	}
	_vias[viaId].addDemand(
			-_trackDemands[layer + 1] * _trackDemands[layer + 1]);
}

void RoutingDB::ripUpVia(const int viaId) {
	reduceViaDemand(viaId);
}

void RoutingDB::restoreVia(const int viaId) {
	addViaDemand(viaId);
}

void RoutingDB::printEdgeInfoRange(const unsigned int start,
		const unsigned int end) {
	for (unsigned i = start; i < end; i++) {
		if (i >= 0 && i < _edges.size()) {
			printEdgeUsage(i);
			printf("\n");
		}
	}
}

void RoutingDB::printViaInfoRange(const unsigned int start,
		const unsigned int end) {
	for (unsigned i = start; i < end; i++) {
		if (i >= 0 && i < _vias.size()) {
			printViaUsage(i);
			printf("\n");
		}
	}
}

