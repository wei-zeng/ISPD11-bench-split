#include "parser.h"
#include <regex>
#include <fstream>

void Parser::ReadAux(const char *auxFile) {
	ReadAux(string(auxFile));
}

void Parser::ReadAux(const string auxFile) {
	const string &fullPathName = auxFile;
	unsigned found = fullPathName.find_last_of("/\\");
	const string path(fullPathName.substr(0, found));

	ifstream inp(fullPathName.c_str());
	if (!inp.is_open())
		ErrorOpenFile(auxFile);

	inp >> tmp;
	inp >> tmp;

	inp >> _nodesFile;
	_nodesFile = path + "/" + _nodesFile;

	inp >> _netsFile;
	_netsFile = path + "/" + _netsFile;

	inp >> _wtsFile;
	_wtsFile = path + "/" + _wtsFile;

	inp >> _plFile;
	_plFile = path + "/" + _plFile;

	inp >> _sclFile;
	_sclFile = path + "/" + _sclFile;

	inp >> _shapesFile;
	_shapesFile = path + "/" + _shapesFile;

	inp >> _routeFile;
	_routeFile = path + "/" + _routeFile;

	cout << _nodesFile << endl;
	cout << _netsFile << endl;
	cout << _plFile << endl;
	cout << _shapesFile << endl;
	cout << _routeFile << endl;

	inp.close();
}

/*void Parser::WriteDEF(const LayoutDR &layout, const RoutingDB_DR &routingDB, const char *defFile) {
	WriteDEF(layout, routingDB, string(defFile));
}

void Parser::WriteDEF(const LayoutDR &layout, const RoutingDB_DR &routingDB, const string defFile) {
	_defwtr.write_def(layout._ldp.def_, routingDB, defFile);
}

void Parser::ReadLEFDEF(LayoutDR &layout, const char *lefFile,
		const char *defFile) {
	ReadLEFDEF(layout, string(lefFile), string(defFile));
}

void Parser::ReadLEFDEF(LayoutDR &layout, const string lefFile,
		const string defFile) {
	layout.readLEFDEF(lefFile, defFile);
}
*/
void Parser::ReadNode(Layout &layout) {
	ifstream inp(_nodesFile);
	if (!inp.is_open())
		ErrorOpenFile(_nodesFile);

	unsigned numNodes = 0;
	unsigned numTerminals = 0;
	while (getline(inp, line)) {
		if (line.empty() || line.at(0) == '#'
				|| line.find("UCLA") != string::npos)
			continue;
		if (line.find("NumNodes") != string::npos) {
			sscanf(line.c_str(), "NumNodes      :  %u", &numNodes);
			continue;
		}
		if (line.find("NumTerminals") != string::npos) {
			sscanf(line.c_str(), "NumTerminals  :  %u", &numTerminals);
			continue;
		}
		layout.addNode(line);
	}
	printf("<I> %-20s : %u\n", "# Nodes", unsigned(layout._nodes.size()));
}
/*
void Parser::ReadNode(LayoutDR &layout) {
	auto &ldp = layout._ldp;
	string objName;
	double width, height, llx, lly;
	string layer;
	int orient;
	for (auto comp : ldp.def_.get_component_umap()) {
		objName = comp.second->name_;
		switch (comp.second->orient_) {
		case 0: // N
		case 2: // S
		case 4: // FN
		case 6: // FS
			width = comp.second->lef_macro_->size_x_;
			height = comp.second->lef_macro_->size_y_;
			break;
		case 1: // W
		case 3: // E
		case 5: // FW
		case 7: // FE
			width = comp.second->lef_macro_->size_y_;
			height = comp.second->lef_macro_->size_x_;
			break;
		}
		llx = comp.second->x_;
		lly = comp.second->y_;
		layer = "metal1";
		orient = comp.second->orient_;
		Node node(objName, width, height, llx, lly, layer, orient);
		layout._nodes.push_back(node);
		layout._nodeUmap.insert(make_pair(objName, layout._numRegularNodes));
		++layout._numRegularNodes;
	}
	for (auto pin : ldp.def_.get_pin_umap()) {
		objName = pin.second->name_;
		width = pin.second->ux_ - pin.second->ly_;
		height = pin.second->uy_ - pin.second->ly_;
		llx = pin.second->lx_;
		lly = pin.second->ly_;
		layer = pin.second->layer_;
		orient = pin.second->orient_;
		Node node(objName, width, height, llx, lly, layer, orient);
		layout._nodes.push_back(node);
		layout._nodeUmap.insert(
				make_pair(objName,
						layout._numRegularNodes + layout._numPTerms));
		++layout._numPTerms;
	}
	printf("<I> %-20s : %u\n", "# Nodes", unsigned(layout._nodes.size()));
}
*/
void Parser::ReadPlace(Layout &layout) {
	ifstream inp(_plFile);
	if (!inp.is_open())
		ErrorOpenFile(_plFile);

	while (getline(inp, line)) {
		if (line.empty() || line.at(0) == '#'
				|| line.find("UCLA") != string::npos)
			continue;
		layout.readNodePlace(line);
	}
}

void Parser::ReadNet(Layout &layout) {
	ifstream inp(_netsFile);
	if (!inp.is_open())
		ErrorOpenFile(_netsFile);

	unsigned numNets = 0;
	unsigned numPins = 0;
	unsigned cntNet = 0;
	while (getline(inp, line)) {
		if (line.empty() || line.at(0) == '#'
				|| line.find("UCLA") != string::npos)
			continue;
		if (line.find("NumNets") != string::npos) {
			sscanf(line.c_str(), "NumNets  :  %u", &numNets);
			layout._nets.resize(numNets);
			continue;
		}
		if (line.find("NumPins") != string::npos) {
			sscanf(line.c_str(), "NumPins  :  %u", &numPins);
			continue;
		}
		if (line.find("NetDegree") != string::npos) {
			numPins = layout._nets[cntNet].init(line);
			for (unsigned i = 0; i < numPins; i++) {
				getline(inp, line);
				layout.addNetPin(line, layout._nets[cntNet]);
			}
			cntNet++;
		}
	}
	printf("<I> %-20s : %u\n", "# Nets", unsigned(layout._nets.size()));
}
/*
void Parser::ReadNet(LayoutDR &layout) {
	auto &ldp = layout._ldp;
	unsigned numNets = ldp.def_.get_net_umap().size();
	unsigned numPins = 0;
	unsigned cntNet = 0;

	layout._nets.resize(numNets);
	for (auto net : ldp.def_.get_net_umap()) {
		numPins += net.second->connections_.size();
		layout._nets[cntNet].setName(net.second->name_);
		for (unsigned i = 0; i < net.second->connections_.size(); i++) {
			layout.addNetPin(net.second->connections_[i], layout._nets[cntNet]);
		}
		cntNet++;
	}
	cout << "# Nets: " << numNets << endl;
	cout << "# Pins: " << numPins << endl;
}
*/
void Parser::ReadShape(Layout &layout) {
	ifstream inp(_shapesFile);
	if (!inp.is_open())
		ErrorOpenFile(_shapesFile);

	unsigned numNonRect = 0;
	unsigned numShapes = 0;

	while (getline(inp, line)) {
		if (line.empty() || line.at(0) == '#'
				|| line.find("UCLA") != string::npos)
			continue;
		if (line.find("NumNonRect") != string::npos) {
			sscanf(line.c_str(), "NumNonRectangularNodes  :  %u", &numNonRect);
			continue;
		}
		if (line.find(":") != string::npos) {
			stringstream ss(line);
			ss >> tmp;
			Node &node = layout._nodes[layout.findNodeId(tmp)];
			node.setIsNonRect();
			ss >> tmp;
			ss >> numShapes;
			for (unsigned i = 0; i < numShapes; i++) {
				getline(inp, line);
				node.addShape(line);
			}
		}
	}
}

void Parser::ReadRouteConfig(Layout &layout) {
	ifstream inp(_routeFile);
	if (!inp.is_open())
		ErrorOpenFile(_routeFile);

	for (unsigned i = 0; i < 4; i++) {
		getline(inp, line);
	}

	getline(inp, line);
	layout.setGrid(line);

	getline(inp, line);
	layout.setVCaps(line);

	getline(inp, line);
	layout.setHCaps(line);

	getline(inp, line);
	layout.setMinWireWidth(line);

	getline(inp, line);
	layout.setMinWireSpacing(line);

	getline(inp, line);
	layout.setViaSpacing(line);

	getline(inp, line);
	layout.setOrig(line);

	getline(inp, line);
	layout.setTileSize(line);

	getline(inp, line);
	layout.setBlkPorosity(line);

	getline(inp, line);
	getline(inp, line);
	unsigned numNiTerms = 0;
	sscanf(line.c_str(), "NumNiTerminals  :  %u", &numNiTerms);

	getline(inp, line);
	for (unsigned i = 0; i < numNiTerms; i++) {
		getline(inp, line);
		char name[256];
		unsigned layer;
		sscanf(line.c_str(), "             %s     %u", name, &layer);
		const string nodeName(name);
		layout._nodes[layout.findNodeId(nodeName)].setLayer(layer - 1);
	}

	getline(inp, line);
	getline(inp, line);
	unsigned numBlocks = 0;
	sscanf(line.c_str(), "NumBlockageNodes  :  %u", &numBlocks);
	printf("<I> %-20s : %u\n", "# Blockage Nodes", numBlocks);

	getline(inp, line);
	for (unsigned i = 0; i < numBlocks; i++) {
		getline(inp, line);
		char name[256];
		sscanf(line.c_str(), "             %s     ", name);
		const string nodeName(name);
		unsigned nodeId = layout.findNodeId(nodeName);
		layout._blkId.push_back(nodeId);  //  record blk node Id into vector
		layout._nodes[nodeId].setBlockedLayers(line);
	}
}
/*
void Parser::ReadRouteConfig(LayoutDR &layout, size_t ub_z) {
	layout.setGrid(ub_z);
	layout.setVCaps();
	layout.setHCaps();
	layout.setMinWireWidth();
	layout.setMinWireSpacing();
	layout.setViaSpacing();
}
*/
// Create files for new designs modified in this program
void Parser::exportDesign(Layout &layout, double netsPercent) {
	exportNodes(layout);
	exportPl(layout);
	exportRoute(layout);
	exportNet(layout, netsPercent);
	exportAux(layout);
}

// Export .node file
void Parser::exportNodes(Layout &layout) {

	char outfile[80];
	sprintf(outfile, "%s.modified.nodes", layout.name().c_str());

	ofstream ofs(outfile);
	ifstream ifs(_nodesFile);

	unsigned numTerminals = 0;
	while (getline(ifs, line)) {
		if (line.find("NumNodes") != string::npos) {
			sscanf(line.c_str(), "NumNodes      :  %u", &_origNumNodes);
			ofs << "NumNodes        :   " << layout._nodes.size() << "\n";
		} else if (line.find("NumTerminals") != string::npos) {
			sscanf(line.c_str(), "NumTerminals :  %u", &numTerminals);
			ofs << "NumTerminals    :   "
					<< (numTerminals + (layout._nodes.size() - _origNumNodes))
					<< "\n";
		} else if (line.find("terminal_NI") != string::npos) {
			break;
		} else {
			ofs << line << "\n";
		}
	}

	for (unsigned nid = _origNumNodes; nid < layout._nodes.size(); nid++) {
		Node &blkNode = layout._nodes[nid];
		ofs << "\t" << blkNode.getName() << "\t" << blkNode.getWidth() << "\t"
				<< blkNode.getHeight() << "\tterminal\n";
	}

	ofs << line << "\n"; /* resume the line we stopped */
	while (getline(ifs, line)) {
		if (line.find("terminal_NI") != string::npos) {
			ofs << line << "\n";
		}
	}

	ifs.close();

	ofs.close();
}

// Export .place file
void Parser::exportPl(Layout &layout) {
	char outfile[80];
	sprintf(outfile, "%s.modified.pl", layout.name().c_str());

	ofstream ofs(outfile);
	ifstream ifs(_plFile);

	while (getline(ifs, line)) {
		if (line.find("FIXED_NI") != string::npos) {
			break;
		}
		ofs << line << "\n";
	}

	for (unsigned nid = _origNumNodes; nid < layout._nodes.size(); nid++) {
		Node &blkNode = layout._nodes[nid];
		ofs << "\t" << blkNode.getName() << "\t" << blkNode.getLLX() << "\t"
				<< blkNode.getLLY() << " : N /FIXED\n";
	}

	ofs << line << "\n";
	while (getline(ifs, line)) {
		if (line.find("FIXED_NI") != string::npos) {
			ofs << line << "\n";
		}
	}
	ifs.close();

	ofs.close();
}

// Export .nets file
void Parser::exportNet(Layout &layout, double netsPercent) {
	char outfile[80];
	sprintf(outfile, "%s.modified.nets", layout.name().c_str());

	ofstream ofs(outfile);

	int numNets = layout._nets.size() * (1 - netsPercent);
	int numPins = 0;

	for (int i = 0; i < numNets; i++) {
		numPins += layout._nets.at(i).numPins();
	}

	ofs << "UCLA nets 1.0\n\n" << "NumNets  :  " << numNets << endl
			<< "NumPins  :  " << numPins << endl << endl;

	for (int i = 0; i < numNets; i++) {
		for (auto line : layout._nets.at(i).infoLines()) {
			ofs << line << "\n";
		}
	}

	ofs.close();
}

// Export .route file
void Parser::exportRoute(Layout &layout) {

	char outfile[80];
	sprintf(outfile, "%s.modified.route", layout.name().c_str());

	ofstream ofs(outfile);
	ifstream ifs(_routeFile);

	unsigned origNumBlks = 0;
	while (getline(ifs, line)) {
		if (line.find("NumBlockageNodes") != string::npos) {
			sscanf(line.c_str(), "NumBlockageNodes  :  %u", &origNumBlks);
			ofs << "NumBlockageNodes    :   "
					<< (origNumBlks + (layout._nodes.size() - _origNumNodes))
					<< "\n";
		} else {
			ofs << line << "\n";
		}
	}

	ifs.close();

	for (unsigned nid = _origNumNodes; nid < layout._nodes.size(); nid++) {
		Node &blkNode = layout._nodes[nid];
		ofs << "\t" << blkNode.getName() << "\t1\t4\n";
	}

	ofs.close();
}

// Export .aux file
void Parser::exportAux(Layout &layout) {
	char outfile[80];
	sprintf(outfile, "%s.modified.aux", layout.name().c_str());

	ofstream ofs(outfile);

	ofs << "RowBasedPlacement :\t" << layout.name() << ".modified.nodes\t"
			<< layout.name() << ".modified.nets\t" << layout.name() << ".wts\t"
			<< layout.name() << ".modified.pl\t" << layout.name() << ".scl\t"
			<< layout.name() << ".shapes\t" << layout.name()
			<< ".modified.route\n";

	ofs.close();
}
