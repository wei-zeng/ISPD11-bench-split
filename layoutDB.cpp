#include "layoutDB.h"

//  .nodes
void Layout::addNode(const string &line) {
	Node node(line);
	_nodes.push_back(node);
	if (node.getName().at(0) == 'o')
		++_numRegularNodes; // movable or terminal
	if (node.getName().at(0) == 'p')
		++_numPTerms; // terminal_NI
}

unsigned Layout::findNodeId(const string &line) const {
	unsigned nodeId;
	if (line.empty()) {
		printf("ERROR %s %s %d\n", __FILE__, __func__, __LINE__);
		printf("Empty node name.\n");
		exit(1);
	}
	if (line.at(0) == 'o') {
		sscanf(line.c_str(), "o%u", &nodeId);
	} else if (line.at(0) == 'p') {
		sscanf(line.c_str(), "p%u", &nodeId);
		nodeId += _numRegularNodes;
	} else {
		printf("ERROR %s %s %d\n", __FILE__, __func__, __LINE__);
		printf("node name : %s\n", line.c_str());
		exit(1);
	}
	return nodeId;
}

//  .nets
void Layout::addNetPin(const string &line, Net &net) {
	istringstream ss(line);

	string tmp;
	ss >> tmp;  //  node name
	int nodeId = findNodeId(tmp);
	Node &node = _nodes[nodeId];
	string pinType;
	ss >> pinType;

	ss >> tmp;
	double offsetX = 0;
	double offsetY = 0;
	ss >> offsetX;
	ss >> offsetY;
	net._netPins.push_back(
			NetPin(node.getPinX(offsetX), node.getPinY(offsetY),
					node.getLayerId()));
	net._netPinTypes.push_back(pinType == "I" ? CI : pinType == "O" ? CO : N);
	net._netPinNodeIDs.push_back(nodeId);
	net._netInfoLines.push_back(line);
	if (pinType == "I") {
		node.addInput();
	} else if (pinType == "O") {
		node.addOutput();
	}
}

//  .nodes
void Layout::readNodePlace(const string &line) {
	istringstream ss(line);
	string name;
	ss >> name;
	unsigned nodeId = findNodeId(name);
	_nodes[nodeId].readPlace(line);
}

//  .route 
void Layout::setGrid(const string &line) {
	sscanf(line.c_str(), "Grid                      :  %u %u %u", &_numTilesX,
			&_numTilesY, &_numLayers);
	printf("<I> %-20s : %u\n", "Grid Dimension X", _numTilesX);
	printf("<I> %-20s : %u\n", "Grid Dimension Y", _numTilesY);
	printf("<I> %-20s : %u\n", "Metal Layout Num", _numLayers);
	_vCaps.resize(_numLayers);
	_hCaps.resize(_numLayers);
	_mww.resize(_numLayers);
	_mws.resize(_numLayers);
	_mvs.resize(_numLayers);
}

void Layout::setVCaps(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	for (unsigned i = 0; i < _numLayers; i++)
		ss >> _vCaps[i];
	printf("<I> %-20s : ", "V Capacity");
	for (unsigned i = 0; i < _numLayers; i++)
		printf("%6.1f", _vCaps[i]);
	printf("\n");
}

void Layout::setHCaps(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	for (unsigned i = 0; i < _numLayers; i++)
		ss >> _hCaps[i];
	printf("<I> %-20s : ", "H Capacity");
	for (unsigned i = 0; i < _numLayers; i++)
		printf("%6.1f", _hCaps[i]);
	printf("\n");
}

void Layout::setMinWireWidth(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	for (unsigned i = 0; i < _numLayers; i++)
		ss >> _mww[i];
	printf("<I> %-20s : ", "Min Wire Width");
	for (unsigned i = 0; i < _numLayers; i++)
		printf("%6.1f", _mww[i]);
	printf("\n");
}

void Layout::setMinWireSpacing(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	for (unsigned i = 0; i < _numLayers; i++)
		ss >> _mws[i];
	printf("<I> %-20s : ", "Min Wire Spacing");
	for (unsigned i = 0; i < _numLayers; i++)
		printf("%6.1f", _mws[i]);
	printf("\n");
}

void Layout::setViaSpacing(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	for (unsigned i = 0; i < _numLayers; i++)
		ss >> _mvs[i];
	printf("<I> %-20s : ", "Min Via Spacing");
	for (unsigned i = 0; i < _numLayers; i++)
		printf("%6.1f", _mvs[i]);
	printf("\n");
}

void Layout::setOrig(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	ss >> _origX;
	ss >> _origY;
}

void Layout::setTileSize(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	ss >> _cellWidth;
	ss >> _cellHeight;
}

void Layout::setBlkPorosity(const string &line) {
	istringstream ss(line);
	string tmp;
	ss >> tmp;
	ss >> tmp;
	ss >> _blkPorosity;
}

void Layout::addNodeBlock(const string &objName, const double width,
		const double height, const double llx, const double lly, const int z) {
	Node nodeBlock(objName, width, height, llx, lly, z, 0);
	_blkId.push_back(_nodes.size());
	_nodes.push_back(nodeBlock);
}

double Layout::getSizeX() const {
    return _cellWidth * _numTilesX;
}

double Layout::getSizeY() const {
    return _cellHeight * _numTilesY;
}
