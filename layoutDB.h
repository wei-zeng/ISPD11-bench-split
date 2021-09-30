#ifndef _LAYOUT_H_
#define _LAYOUT_H_

#include "baseDB.h"
#include "nodes.h"
#include "nets.h"
#include "vias.h"
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Layout {
	friend class Parser;
	friend class RoutingDB;
public:
	Layout() :
			_numTilesX(0), _numTilesY(0), _numLayers(0), _numRegularNodes(0), _numPTerms(
					0) {
	}

	void initDesignName(const string &designName) {
		_designName = designName;
	}
	void initDesignName(const char *designName) {
		_designName = string(designName);
	}
	string name() const {
		return _designName;
	}
	//  .nodes
	void addNode(const string &line);
	unsigned findNodeId(const string &line) const;

	//  .nets
	void addNetPin(const string &line, Net &net);

	//  .pl
	void readNodePlace(const string &line);

	//  .shapes

	//  .route
	void setGrid(const string &line);
	void setVCaps(const string &line);
	void setHCaps(const string &line);
	void setMinWireWidth(const string &line);
	void setMinWireSpacing(const string &line);
	void setViaSpacing(const string &line);
	void setOrig(const string &line);
	void setTileSize(const string &line);
	void setBlkPorosity(const string &line);
	double getCellWidth() {return _cellWidth;}
	double getCellHeight() {return _cellHeight;}
    double getSizeX() const;
    double getSizeY() const;

	//  modifier
	void addNodeBlock(const string &objName, const double width,
			const double height, const double llx, const double lly,
			const int z);

protected:

	string _designName;
	unsigned _numTilesX;
	unsigned _numTilesY;
	unsigned _numLayers;
	//  .nodes file and .pl file
	vector<Node> _nodes;
	unsigned _numRegularNodes;
	unsigned _numPTerms;

	//  .nets file
	vector<Net> _nets;

	//  .route file
	vector<double> _vCaps;
	vector<double> _hCaps;
	vector<double> _mww;   //  minimum wire width
	vector<double> _mws;   //  minimum wire spacing
	vector<double> _mvs;   //  minimum via spacing
	double _origX, _origY;
	double _cellWidth, _cellHeight;
	double _blkPorosity;

	//  blockage
	vector<unsigned> _blkId;

	//  global nets
	Gcell getGcell(const Point3D &pt) const {
		return Gcell(static_cast<int>((pt._x - _origX) / _cellWidth),
				static_cast<int>((pt._y - _origY) / _cellHeight), pt._z);
	}
};

#endif
