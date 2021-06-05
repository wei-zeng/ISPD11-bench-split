#ifndef _NETS_H_
#define _NETS_H_

#include "baseDB.h"
#include "vias.h"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>

typedef Point3D NetPin;
typedef Gcell GnetPin;

class Net {
	friend class Layout;
	//friend class LayoutDR;
	friend class RoutingDB;
	//friend class RoutingDB_DR;
public:
	Net() {
	}
	int init(const string &line) {
		_netInfoLines.push_back(line);
		istringstream ss(line);
		string tmp;
		int numPins;
		ss >> tmp;
		ss >> tmp;
		ss >> numPins;
		_netPins.reserve(numPins);
		ss >> _name;
		return numPins;
	}

	string getName() const {
		return _name;
	}
	void setName(const string &name) {
		_name = name;
	}
	int numPins() const {
		return _netPins.size();
	}
	vector<string> infoLines() const {
		return _netInfoLines;
	}

protected:
	string _name;
	vector<NetPin> _netPins;
	vector<PinType> _netPinTypes;
	vector<int> _netPinNodeIDs;
	vector<string> _netInfoLines;
};

class Gnet {
	friend class RoutingDB;
	//friend class RoutingDB_DR;
	//friend class MyDefWriter;
public:
	Gnet() :
			_isGlobal(false) {
	}
	bool isGlobal() const {
		return _isGlobal;
	}
	bool isRealNet() const {
		return _isRealNet;
	}
	int parent(const int curr) const {
		return _gWires[curr]._pWireId;
	}
	int parent2d(const int curr) const {
		return _wire2d[curr]._pWireId;
	}
	string getName() const {
		return _name;
	}
	bool operator==(const Gnet &other) const {
		return (_gWires == other._gWires);
	}
	void printGnet() const {
		cout << _name << endl;
		for (auto gWireItr1 = _gWires.begin(); gWireItr1 < _gWires.end();
				++gWireItr1) {
			cout << gWireItr1 - _gWires.begin() << " ";
			(*gWireItr1).printGlobalWire();
			cout << flush;
		}
	}
	void removeEdge(int e) {
		_removedEdges.push_back(e);
	}
	void removeVia(int v) {
		_removedVias.push_back(v);
	}
	void addEdge(int e) {
		_addedEdges.push_back(e);
	}
	void addVia(int v) {
		_addedVias.push_back(v);
	}
	vector<int> getRemovedEdges() const {
		return _removedEdges;
	}
	vector<int> getRemovedVias() const {
		return _removedVias;
	}
	bool findCleanEdge(size_t edgeId) const {
		return _cleanEdges.find(edgeId) != _cleanEdges.end();
	}
	bool findCleanVia(size_t viaId) const {
		return _cleanVias.find(viaId) != _cleanVias.end();
	}
	void setCleanEdge(size_t edgeId, const Edge edge) {
		_cleanEdges.insert(make_pair(edgeId, edge));
	}
	void setCleanVia(size_t viaId, const Via via) {
		_cleanVias.insert(make_pair(viaId, via));
	}
	void resetCleanEdge() {
		_cleanEdges.clear();
	}
	void resetCleanVia() {
		_cleanVias.clear();
	}
	bool findRUEdge(size_t edgeId) const {
		return _ruEdges.find(edgeId) != _ruEdges.end();
	}
	bool findRUVia(size_t viaId) const {
		return _ruVias.find(viaId) != _ruVias.end();
	}
	void resetRUEdge() {
		_ruEdges.clear();
	}
	void resetRUVia() {
		return _ruVias.clear();
	}
	void setRUEdge(size_t edgeId, const Edge edge) {
		_ruEdges.insert(make_pair(edgeId, edge));
	}
	void setRUVia(size_t viaId, const Via via) {
		_ruVias.insert(make_pair(viaId, via));
	}
private:
	//  add point pair
	void addPtPair(const Gcell &gp1, const Gcell &gp2) {
		_gpPairArray.push_back(make_pair(gp1, gp2));
	}
	bool _isGlobal;
	bool _isRealNet;
	size_t _root;
	vector<GnetPin> _gPins;
	vector<GnetPin> _gPins2D;
	string _name;
	Net _net;
	//  read from global wires
	vector<pair<Gcell, Gcell> > _gpPairArray;
	vector<GlobalWire> _gWires;
	vector<PinType> _pinTypes;
	vector<Node> _nodes;

	//  convert from global wires to 2D wires
	vector<Wire2d> _wire2d;
	vector<int> _removedEdges, _addedEdges;
	vector<int> _removedVias, _addedVias;
	map<size_t, Edge> _cleanEdges;
	map<size_t, Via> _cleanVias;
	map<size_t, Edge> _ruEdges;
	map<size_t, Via> _ruVias;
	set<size_t> _occupiedEdges;
	set<size_t> _occupiedVias;
};

#endif
