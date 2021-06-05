#ifndef _NODES_H_
#define _NODES_H_

#include "baseDB.h"
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

class Node {
	friend class RoutingDB;
	//friend class RoutingDB_DR;
public:
	Node() :
			_isTerm(false), _isNonRect(false), _layerId(0), _numInputs(0), _numOutputs(
					0), _orient(0) {
	}
	Node(const string &line) :
			_isTerm(false), _isNonRect(false), _layerId(0), _numInputs(0), _numOutputs(
					0), _orient(0) {
		init(line);
	}
	Node(const string &objName, const double width, const double height,
			const double llx, const double lly, const string &layer, const int orient) :
			_name(objName), _width(width), _height(height), _llx(llx), _lly(
					lly), _isTerm(false), _isNonRect(false), _numInputs(0), _numOutputs(
					0), _orient(orient) {
		int z = stoi(layer.substr(5)) - 1;
		_layerId = z;
	}
	Node(const string &objName, const double width, const double height,
			const double llx, const double lly, const int z, const int orient) :
			_name(objName), _width(width), _height(height), _layerId(z), _llx(
					llx), _lly(lly), _isTerm(false), _isNonRect(false), _numInputs(
					0), _numOutputs(0), _orient(orient) {
		_blockedLayers.push_back(z);
	}

	void init(const string &line) {
		istringstream ss(line);
		string tmp;
		ss >> _name;
		ss >> _width;
		ss >> _height;
		ss >> tmp;
		if (tmp.compare("terminal") == 0 || tmp.compare("terminal_NI") == 0) {
			_isTerm = true;
		}
	}

	void readPlace(const string &line) {
		istringstream ss(line);
		string name;
		ss >> name;
		ss >> _llx;
		ss >> _lly;
		_layerId = 0;
	}

	void addShape(const string &line) {
		istringstream ss(line);
		string tmp;
		ss >> tmp;
		double llx, lly, width, height;
		ss >> llx;
		ss >> lly;
		ss >> width;
		ss >> height;
		_shapes.push_back(Shape(llx, lly, width, height));
	}

	void setLayer(const unsigned layerId) {
		_layerId = layerId;
	}
	void setIsNonRect() {
		_isNonRect = true;
	}

	void setBlockedLayers(const string &line) {
		istringstream ss(line);
		string tmp;
		ss >> tmp;
		unsigned num, layer;
		ss >> num;
		for (unsigned i = 0; i < num; i++) {
			ss >> layer;
			_blockedLayers.push_back(layer - 1);
		}
	}

	void addInput() {
		_numInputs++;
	}
	void addOutput() {
		_numOutputs++;
	}

	string getName() const {
		return _name;
	}
	double getPinX(double offsetX) const {
		return double(ceil(_llx + _width / 2 + offsetX));
	}
	double getPinY(double offsetY) const {
		return double(ceil(_lly + _height / 2 + offsetY));
	}
	unsigned getLayerId() const {
		return _layerId;
	}
	bool isNonRect() const {
		return _isNonRect;
	}
	double getLLX() const {
		return _llx;
	}
	double getLLY() const {
		return _lly;
	}
	double getWidth() const {
		return _width;
	}
	double getHeight() const {
		return _height;
	}

	void printNode() const {
		if (_blockedLayers.empty())
			cout << "<D> " << _name << " Regular Blockage" << endl;
		else {
			if (_isNonRect)
				cout << "<D> " << _name << " Non-rect Blockage" << endl;
			else
				cout << "<D> " << _name << " Rect Blockage" << endl;
		}
	}

private:
	string _name;
	double _width, _height;
	bool _isTerm;
	double _llx, _lly;
	unsigned _layerId;
	//  non-rect
	bool _isNonRect;
	int _numInputs, _numOutputs;
	int _orient;
public:
	vector<Shape> _shapes;
	//  block
	vector<int> _blockedLayers;
};

#endif
