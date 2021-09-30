#ifndef _GRDB_H_
#define _GRDB_H_

#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <random>

#include "baseDB.h"
#include "layoutDB.h"
#include "graph.h"

class Vpin {
public:
	int vPinID;
	int gnetID;
	int matchingVpinIdx = -1;
	int xCoord, yCoord, zCoord;
	int wlToL1;
	bool isPIO = false;
	bool hasOriginalRoot = false;
	std::unordered_set<Gcell, HashGcell3d> gnetVertices; // endpoints of gWires
	std::unordered_set<Gcell, HashGcell3d> gnetGrids; // endpoints of gWires and grids in between
    std::vector<std::tuple<int, int, int, int, int, int>> routes;
    std::vector<Point3D> pins;
	std::vector<PinType> pinTypes; // 0=N, 1=CI, 2=CO, 3=PI, 4=PO
	std::vector<int> pinNodeIDs;
	std::vector<double> cellAreas;
	int vpGroupID; // for grouping vpins that share the same node
};

struct ObfusOption {
	int gnetID;
	int d0, d1, axis;
	int margin0, margin1, margin2;
	Gnet newGnet;
	bool operator>(const ObfusOption &other) const;
	std::string toString() const;
};

struct VpinPair {
	int gnetID;
	bool obfuscated = false;
	double gainUpperBound = 2.0;
	double gainLowerBound = -1.0;
};

class RoutingDB {

	//
	//  Layout Information
	//
public:
	RoutingDB(const string &designName, const Layout &layout) :
			_designName(designName), _numTilesX(layout._numTilesX), _numTilesY(
					layout._numTilesY), _numLayers(layout._numLayers), _numGlobalNets(
					0), _tileWidth(static_cast<int>(layout._cellWidth)), _tileHeight(
					static_cast<int>(layout._cellHeight)), _dirty(false) {
	}
	bool operator==(const RoutingDB &other) {
		return _gnets == other._gnets && _edges == other._edges
				&& _vias == other._vias;
	}

private:
	string _designName;
	int _numTilesX, _numTilesY, _numLayers;
	vector<RoutingDir> _dirLayers;
	vector<int> _trackDemands;
	vector<int> _cntEdgeLayers;
	vector<int> _capOnLayer;
	vector<int> _upViaSize;
	vector<int> _dnViaSize;
	int _tileWidth, _tileHeight;
	//
	//  Edge Related Objects
	//
public:
	void initEdges(const Layout &layout);
	void initGlobalEdgeProfile(const Layout &layout);
	void printEdgeInfoRange(const unsigned int start, const unsigned int end);
	void printViaInfoRange(const unsigned int start, const unsigned int end);
private:
	void refreshNodeBlkInfo(const Layout &layout, const Node &node);
	void refreshRectBlkInfo(const Layout &layout, const Shape &shape,
			const int layerId);
	void updateEdgeFreeList(const Layout &layout, const int edgeId,
			const Shape &shape);
	void updateSingleEdgeFreeList(Edge &edge, const double blkStart,
			const double blkEnd);

	// int findEdge( const int x, const int y, const int z ) const;
	RoutingDir getRoutingDir(const int layer) const;
	// Gcell getEdgeLoc(const int edgeId) const;
	int getEdgeLayer(const int edgeId) const;
	// int getEdgeDemand( const int edgeId ) const;
	int getEdgeOf(const int edgeId) const;
	void addEdgeDemand(const int edgeId);
	void reduceEdgeDemand(const int edgeId);
	void ripUpEdge(const int edgeId);
	void restoreEdge(const int edgeId);

public:
	int findEdge(const int x, const int y, const int z) const;
	Gcell getEdgeLoc(const int edgeId) const;
	vector<Edge> getEdges() const {
		return _edges;
	}
	int getEdgeDemand(const int edgeId) const;
	int getEdgeBlk(const int edgeId) const;
	int getEdgeCap(const int edgeId) const;
	void clearEdgeDemands();
private:

	void checkEdgeFreeList(const int edgeId);
	//  gWire related methods
	bool isUpViaGlobalWire(const Gnet &gnet, const GlobalWire &gWire) const;
	bool isDnViaGlobalWire(const Gnet &gnet, const GlobalWire &gWire) const;
	void bruteForceFindAllChildWires(const Gnet &gnet, const int wireId,
			vector<int> &cache);

	inline void rmEdgeDemand(const int edgeId);
	vector<Edge> _edges;
	void printEdgeUsage(const size_t edgeId);
	void printCompleteEdgeInfo(const int edgeId);

	//
	//  Via Related Objects
	//
public:
	void initVias(const Layout &layout);
	vector<Via> getVias() const {
		return _vias;
	}
private:
	long findVia(const int x, const int y, const int z) const;
	void addViaDemand(const int viaId);
	void reduceViaDemand(const int viaId);
	void ripUpVia(const int viaId);
	void restoreVia(const int viaId);

	Gcell getViaLoc(const int viaId) const;
	int getViaLayer(const int viaId) const;
	int getViaOf(const int viaId) const;

	vector<Via> _vias;
	void printViaUsage(const size_t viaId);

	//
	//  Read Global Routes
	//
public:
	void initGlobalNets(Layout &layout);
	void readGlobalWires(const Layout &layout, const char *wiresFile);
	void readGlobalWires(const Layout &layout, const string &wiresFile);
	size_t writeGlobalWires(const Layout &layout, const char *wiresFile) const;
	size_t writeGlobalWires(const Layout &layout, const char *wiresFile, const std::unordered_set<size_t> &excludedGnetIds) const;
	void writeGlobalWiresVpin(const Layout &layout, const char *wiresFile, const Vpin &vp, size_t &netCount) const;
	void writeNets(const Layout &layout, const char *netsFile, const vector<Vpin> &vpins, bool brokenNetsOnly) const;
	void writeKey(const Layout &layout, const char *keyFile, const vector<Vpin> &vpins) const;
	void initRoutingTreeForAllGnets();
	bool checkAllRoutingTrees();
	int getViaDemand(const int viaId) const;
	int getViaCap(const int viaId) const;
	vector<Gnet>& getGnets() {
		return _gnets;
	}
private:
	void initRoutingTree(Gnet &gnet);

	//
	//  Routing Tree objects
	//
	vector<vector<PtInfo> > _treeTable;
	// Structure of PtInfo is in baseDB.h: (int x, int y, int z, int indexId, int distPrev, int distNext, bool hasDnVia, bool hasUpVia)
	vector<PtInfo> tmp;

	vector<unordered_map<Gcell, int, HashGcell2d> > _uMapCell2LocInTable;
	vector<Gcell> frontier;
	vector<Gnet> _gnets;
	int _numGlobalNets;
	/* used for removing antenna wires */
	vector<int> _offset;
	vector<int> _numChildren;
	vector<bool> _validWires;
	vector<GlobalWire> _cleanWires;
    bool _dirty;
    bool checkChildren(Gnet &gnet, vector<int> &_numChildren,
			vector<bool> &_validWires, double &wlReduced, int splitLayer,
			double &wlReducedAboveSL,
			unordered_set<Gcell, HashGcell3d> &remainingGnetVertices,
			unordered_set<Gcell, HashGcell3d> &remainingGnetGrids, int curr,
			const Vpin &vp);

public:
	void updateEdgeDemands();
    std::tuple<vector<Vpin>, vector<int>, vector<int>, double, double> ripUpPin(
			const Layout &layout, Vpin &vp, const int splitLayer);
	std::tuple<double, vector<int>, vector<int>, int,
			std::unordered_set<Gcell, HashGcell3d>> rerouteNet(
			const Layout &layout, Vpin &vp,
			const std::unordered_set<Gcell, HashGcell3d> &bannedPts,
			const int zlLim, const int zuLim, const int maxMargin,
			const bool singleMargin, const bool writeToDB,
			const bool updateRoot);
	std::tuple<double, vector<int>, vector<int>, int> rerouteNet(
			const Layout &layout, const Vpin &vp1, const Vpin &vp2,
			const int zlLim, const int zuLim, const int maxMargin,
			const bool singgleMargin, const bool writeToDB);
	void restoreNet(const int gnetID, const Gnet gnet,
			const std::map<size_t, Edge> &cleanEdges,
			const std::map<size_t, Via> &cleanVias, bool restoreCongestion);
	int getTotalWL() const;
	int getWlByNet(int gnetID) const;
	int getTotalWLByDemand() const;
    bool isDirty() const {
		return _dirty;
	}
    void setDirty() {
        _dirty = true;
    }
	vector<Vpin> getVpins(const Layout &layout, const int splitLayer,
			const bool twoCutNetsOnly, bool includeFloatingVpins, bool excludeNI) const;
private:
	void createBlockage(const int xStart, const int yStart, const int xEnd,
			const int yEnd, const int z, Layout &layout);

// A-star:
public:
	Graph build_graph(const size_t gnetId,
			const std::unordered_set<Gcell, HashGcell3d> &bannedPt,
			const int xlLim, const int xuLim, const int ylLim, const int yuLim,
			const int zlLim, const int zuLim) const;
	unordered_set<Vertex> build_goals(const Vpin &vp, const int xlLim,
			const int xuLim, const int ylLim, const int yuLim,
			const int zlLim) const;

	std::pair<double, vector<Vertex>> a_star(Graph &graph3D, const Vertex &s,
			const unordered_set<Vertex> &goals) const;
	std::pair<ObfusOption, bool> randomTraverse(const Layout &layout, int layer,
			vector<Vpin> &vpins, size_t best_vpin, size_t maxTrials,
			double stdevX, double stdevY, std::mt19937 &rng);
    bool applyChange(const Layout &layout, int layer,
			vector<Vpin> &vpins, size_t best_vpin, const ObfusOption &opt);
};

bool congestionEqual(const RoutingDB &db1, const RoutingDB &db2);

#endif
