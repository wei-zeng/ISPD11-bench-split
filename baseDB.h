#ifndef _BASEDB_H_
#define _BASEDB_H_

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cfloat>
#include<vector>
#include <list>
#include<algorithm>

#define errorExit printf( "ERROR %s %s %d\n", __FILE__, __func__, __LINE__ )

using namespace std;

enum RoutingDir {
	V, H, NA
};

namespace WireDir {
enum Dir {
	NegativeX, PositiveX, NegativeY, PositiveY, Undef
};
}

enum PinType {
	N, CI, CO, PI, PO
};

struct Point2D {
	int _x;
	int _y;
	Point2D(int x, int y) :
			_x(x), _y(y) {
	}
};

class Point3D {
public:
	Point3D() {
	}
	Point3D(double x, double y, int z) :
			_x(x), _y(y), _z(z), _lx(x), _ux(x), _ly(y), _uy(y) {
	}
	/*Point3D(double lx, double ly, double ux, double uy, int z, Polygon shape) :
			_x((lx + ux) / 2), _y((ly + uy) / 2), _z(z), _lx(lx), _ux(ux), _ly(
					ly), _uy(uy), _shape(shape) {
	}*/
	double _x, _y;
	double _lx, _ly, _ux, _uy; // for bbox of pins
	//Polygon _shape; // for pins
	int _z;

	double operator-(const Point3D &other) const {
		return abs(_x - other._x) + abs(_y - other._y);
	}
};

class Gcell {
public:
	Gcell() = default;
	Gcell(int x, int y, int z) :
			_x(x), _y(y), _z(z), _lx(x), _ux(x), _ly(y), _uy(y) {
	}
	/*Gcell(int lx, int ly, int ux, int uy, int z, Polygon shape) :
			_x((lx + ux) / 2), _y((ly + uy) / 2), _z(z), _lx(lx), _ux(ux), _ly(
					ly), _uy(uy), _shape(shape) {
	}*/
	int _x, _y, _z;
	int _lx, _ly, _ux, _uy; // for bbox of gpins
	//Polygon _shape; // for gpins
	bool operator==(const Gcell &other) const {
		return _x == other._x && _y == other._y && _z == other._z;
	}
	bool operator<(const Gcell &other) const {
		return _z != other._z ? _z < other._z :
				_y != other._y ? _y < other._y : _x < other._x;
	}
	int operator-(const Gcell &other) const {
		return abs(_x - other._x) + abs(_y - other._y) + abs(_z - other._z);
	}
	bool overlapWith(const Gcell &other) const {
		return (_lx <= other._ux && _ux >= other._lx && _ly <= other._uy
				&& _uy >= other._ly);
	}
	int viaWeights[9] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	int viaWeightedDist(const Gcell &other) const {
		int dist = abs(_x - other._x) + abs(_y - other._y);
		for (int z = min(_z, other._z); z < max(_z, other._z); z++)
			dist += viaWeights[z];
		return dist;
	}
};

class Shape {
public:
	Shape(double llx, double lly, double width, double height) :
			_llx(llx), _lly(lly), _width(width), _height(height) {
	}
	double _llx, _lly, _width, _height;
};

struct GlobalWire {
	GlobalWire() {
	}
	GlobalWire(const int x, const int y, const int z, const int pWireId,
			const int realPinId) :
			_x(x), _y(y), _z(z), _pWireId(pWireId), _realPinId(realPinId) {
	}

	bool operator<(const GlobalWire &other) const {
		return _z != other._z ? _z < other._z :
				_y != other._y ? _y < other._y : _x < other._x;
	}
	bool operator==(const GlobalWire &other) const {
		return _x == other._x && _y == other._y && _z == other._z
				&& _pWireId == other._pWireId && _realPinId == other._realPinId;
	}
	int _x, _y, _z;
	int _pWireId;
	int _realPinId;
	virtual void printGlobalWire() const {
		printf("gWire : %d, %d, %d, %d, %d\n", _x, _y, _z, _pWireId,
				_realPinId);
	}
};

struct GcellUsageInfo {
	GcellUsageInfo() :
			_negativeX(false), _positiveX(false), _negativeY(false), _positiveY(
					false), _isRealPin(false), _wireIndexInOutputWireVector(-1), _score(
					0) {
	}
	void setupstreamGcellUsageInfo(WireDir::Dir dir) {
		if (dir == WireDir::NegativeX)
			_positiveX = true;
		if (dir == WireDir::NegativeY)
			_positiveY = true;
		if (dir == WireDir::PositiveX)
			_negativeX = true;
		if (dir == WireDir::PositiveY)
			_negativeY = true;
	}

	bool connected() const {
		return _negativeX || _positiveX || _negativeY || _positiveY;
	}
	bool isRealPin() const {
		return _isRealPin;
	}
	bool _negativeX, _positiveX, _negativeY, _positiveY, _isRealPin;
	int _wireIndexInOutputWireVector;
	float _score;
};

class Wire2d: public GlobalWire {
	friend class RoutingDB;
public:
	Wire2d(const int x, const int y, const int pWireId, const int realPinId) :
			GlobalWire(x, y, 0, pWireId, realPinId), _len(-1) {
		_realPinArray.assign(16, false);
	}
	void printGlobalWire() const {
		printf("wire2d : x=%d, y=%d, parent=%d, pinId=%d, len=%d, dir=%s\n", _x,
				_y, _pWireId, _realPinId, _len,
				_dir == WireDir::NegativeX ? "negX" :
				_dir == WireDir::PositiveX ? "posX " :
				_dir == WireDir::NegativeY ? "negY" :
				_dir == WireDir::PositiveY ? "posY " : "undef");
	}

private:
	WireDir::Dir _dir;
	int _len;    //  the length of the wire2d (i.e. edge number)
	vector<int> _layers; //  has as size of _length; initialized as -1; from source to parent
	vector<GcellUsageInfo> _gcellsInWire;    //  len * numLayers
	vector<int> _childWires;
	vector<bool> _realPinArray;
};

class Edge {
public:
	Edge() :
			_cap(0), _blk(0), _demand(0) {
	}

//  write
	void setCap(int cap) {
		_cap = cap;
	}
	void setBlk(int blk) {
		_blk = blk;
	}
	void setLoc(int x, int y, int z) {
		_x = x;
		_y = y;
		_z = z;
	}
	void addDemand(int changeDemand) {
		_demand += changeDemand;
	}
	void clearDemand() {
		_demand = 0;
	}
	void ripUp() {
		_rippedUp = true;
	}
	void restore() {
		_rippedUp = false;
	}
	void markOriginal() {
		_original = true;
	}
	void demarkOriginal() {
		_original = false;
	}

// read
	int cap() const {
		return _cap;
	}
	int blk() const {
		return _blk;
	}
	int demand() const {
		return _demand;
	}

	int of() const {
		return max(0, _demand + _blk - _cap);
	}
	Gcell getLoc() const {
		return Gcell(_x, _y, _z);
	}
	bool isRippedUp() const {
		return _rippedUp;
	}
	bool isOriginal() const {
		return _original;
	}

	void initFreeList(const double start, const double end) {
		_freeList.insert(_freeList.begin(), make_pair(start, end));
	}

	bool operator==(const Edge &other) const {
		return _blk == other._blk && _cap == other._cap
				&& _demand == other._demand;
	}

	double getFreeLength() const {
		double available = 0;
		for (auto seg : _freeList)
			available += (seg.second - seg.first);
		return available;
	}

//  debug
	void printEdgeUsage() const {
		printf("cap : %-10dblk : %-10ddemand : %-10d\n", _cap, _blk, _demand);
	}

	list<pair<double, double> > _freeList;

private:

//  global routing
	int _cap;
	int _blk;
	int _demand;
	bool _rippedUp;
	bool _original;
	int _x, _y, _z; // DR only
};

struct PtInfo {
	PtInfo(int x, int y, int z, int indexId, int distPrev, int distNext,
			bool hasDnVia, bool hasUpVia) :
			_x(x), _y(y), _z(z), _visited(false), _indexId(indexId), _distPrev(
					distPrev), _distNext(distNext), _distPrevNP(-1), _distNextNP(
					-1), _hasDnVia(hasDnVia), _hasUpVia(hasUpVia) {
	}
	PtInfo(int x, int y, int z, int indexId, int distPrev, int distNext,
			int distPrevNP, int distNextNP, bool hasDnVia, bool hasUpVia) :
			_x(x), _y(y), _z(z), _visited(false), _indexId(indexId), _distPrev(
					distPrev), _distNext(distNext), _distPrevNP(distPrevNP), _distNextNP(
					distNextNP), _hasDnVia(hasDnVia), _hasUpVia(hasUpVia) {
	}

	bool operator==(const PtInfo &other) const {
		return _x == other._x && _y == other._y && _z == other._z;
	}
	bool operator!=(const PtInfo &other) const {
		return _x != other._x || _y != other._y || _z != other._z;
	}

	int operator-(const PtInfo &other) const {
		return abs(_x - other._x) + abs(_y - other._y) + abs(_z - other._z);
	}

	int _x, _y, _z, _indexId;
	bool _visited, _hasDnVia, _hasUpVia;
	int _distPrev, _distNext;
	int _distPrevNP, _distNextNP; // in non-preferred direction

	Gcell getLoc() const {
		return Gcell(_x, _y, _z);
	}
	void printPtInfo() {
		printf("%4d,%4d,%4d,%4d,%4d,%4d,%s,%s\n", _x, _y, _z, _indexId,
				_distPrev, _distNext, _hasDnVia ? "TRUE" : "FALSE",
				_hasUpVia ? "TRUE" : "FALSE");
	}
};

struct comparePtInfoOnHorLayer {
	bool operator()(const PtInfo &pt1, const PtInfo &pt2) {
		return pt1._y != pt2._y ? pt1._y < pt2._y : pt1._x < pt2._x;
	}
};

struct comparePtInfoOnVerLayer {
	bool operator()(const PtInfo &pt1, const PtInfo &pt2) {
		return pt1._x != pt2._x ? pt1._x < pt2._x : pt1._y < pt2._y;
	}
};

struct HashGcell2d {
	size_t operator()(const Gcell &pin) const {
		return (pin._x << 16) | pin._y;
	}
};

struct HashGcell3d {
	size_t operator()(const Gcell &pin) const {
		return (pin._z << 28) | (pin._y << 14) | pin._x;
	}
};

#endif
