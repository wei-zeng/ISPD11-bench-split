#ifndef _PARSER_H_
#define _PARSER_H_

#include "grDB.h"
//#include "grDB_DR.h"
#include "graph.h"
//#include "DefWriter.h"
#include <cstring>
#include <iostream>

inline void ErrorOpenFile(const string &fileName) {
	cout << "Cannot open file : " << fileName << endl;
	exit(1);
}

class Parser {
public:
	Parser() {}//:
			//_defwtr(MyDefWriter::get_instance()) {
	void ReadAux(const string auxFile);
	void ReadAux(const char *auxFile);
	//void ReadLEFDEF(LayoutDR &layout, const string lefFile,
	//		const string defFile);
	//void ReadLEFDEF(LayoutDR &layout, const char *lefFile, const char *defFile);
	void ReadNode(Layout &layout);
	//void ReadNode(LayoutDR &layout);
	void ReadPlace(Layout &layout);
	void ReadNet(Layout &layout);
	//void ReadNet(LayoutDR &layout);
	void ReadShape(Layout &layout);
	//void ReadShape(LayoutDR &layout);
	void ReadRouteConfig(Layout &layout);
	//void ReadRouteConfig(LayoutDR &layout, size_t ub_z);
	void exportDesign(Layout &layout, double netsPercent);
	//void WriteDEF(const LayoutDR &layout, const RoutingDB_DR &routingDB,
	//		const string defFile);
	//void WriteDEF(const LayoutDR &layout, const RoutingDB_DR &routingDB,
	//		const char *defFile);

private:
	void exportNodes(Layout &layout);
	void exportPl(Layout &layout);
	void exportRoute(Layout &layout);
	void exportNet(Layout &layout, double netsPercent);
	void exportAux(Layout &layout);

	string _nodesFile;
	string _netsFile;
	string _wtsFile;
	string _plFile;
	string _sclFile;
	string _shapesFile;
	string _routeFile;
	//MyDefWriter &_defwtr;

	string line, tmp;

	unsigned _origNumNodes;
};

#endif
