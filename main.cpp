/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  Main file
 *
 *        Version:  1.0
 *        Created:  11/01/2014 16:56:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Daohang Shi, Jackson Melchert, Wei Zeng
 *
 * =====================================================================================
 */
#include <cstdio>
#include <sys/param.h>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <string>
#include "parser.h"
#include <map>
#include <algorithm>
#include <numeric>

int main(int argc, char **argv) {
	char *designName = nullptr;
	char *auxFile = nullptr;
	char *rtFile = nullptr;
	char *outputRtFile = nullptr;
	char *outputNetsFile = nullptr;
	char *outputKeyFile = nullptr;
	int i = 1;
	int layer = 0;
    bool brokenNetsOnly = false;
    bool twoCutNetsOnly = false;
	bool floatingNets = false; // include floating nets (not connected to any pin)
	bool excludeNI = false; // exclude FIXED_NI terminals ("p" nodes in .nodes files) in vpins
    bool obfuscate = false; // enable random obfuscation 
	// Read input arguments
	while (i < argc) {
		if (!strcmp(argv[i], "-design")) {
			i++;
			assert(i < argc);
			designName = argv[i++];
        } else if (!strcmp(argv[i], "-layer")) {
			i++;
			assert(i < argc);
			layer = stoi(argv[i++]);
		} else if (!strcmp(argv[i], "-auxFile")) {
			i++;
			assert(i < argc);
			auxFile = argv[i++];
		} else if (!strcmp(argv[i], "-rtFile")) {
			i++;
			assert(i < argc);
			rtFile = argv[i++];
		} else if (!strcmp(argv[i], "-outputRt")) {
			i++;
			assert(i < argc);
			outputRtFile = argv[i++];
		} else if (!strcmp(argv[i], "-outputNets")) {
			i++;
			assert(i < argc);
			outputNetsFile = argv[i++];
		} else if (!strcmp(argv[i], "-outputKey")) {
			i++;
			assert(i < argc);
			outputKeyFile = argv[i++];
        } else if (!strcmp(argv[i], "-brokenNetsOnly")) {
            i++;
            brokenNetsOnly = true;
        } else if (!strcmp(argv[i], "-twoCutNetsOnly")) {
            i++;
            brokenNetsOnly = true;
            twoCutNetsOnly = true;
        } else if (!strcmp(argv[i], "-floatingNets")) {
            i++;
            floatingNets = true;
        } else if (!strcmp(argv[i], "-excludeNI")) {
            i++;
            excludeNI = true;
        } else if (!strcmp(argv[i], "-obfuscate")) {
            i++;
            obfuscate = true;
		} else {
			cerr << "Usage : " << argv[0] << " -design <design_name> -layer <split_layer> -auxFile <aux_file> -rtFile <rt_file> -outputRt <output_rt_file> -outputNets <output_nets_file> -outputKey <output_key_file> [-brokenNetsOnly] [-twoCutNetsOnly] [-floatingNets] [-excludeNI] [-obfuscate]" << endl;
			exit(1);
		}
	}

	if (!designName || !auxFile || !rtFile || layer == 0 ||
            !outputRtFile || !outputNetsFile || !outputKeyFile) {
		cerr << "Usage : " << argv[0] << " -design <design_name> -layer <split_layer> -auxFile <aux_file> -rtFile <rt_file> -outputRt <output_rt_file> -outputNets <output_nets_file> -outputKey <output_key_file> [-brokenNetsOnly] [-twoCutNetsOnly] [-floatingNets] [-excludeNI] [-obfuscate]" << endl;
		exit(1);
	}

	Parser parser;
	Layout layout;
	layout.initDesignName(designName);

	// Read design
	if (auxFile) {
		parser.ReadAux(auxFile);
		parser.ReadNode(layout);
		parser.ReadPlace(layout);
		parser.ReadShape(layout);
		parser.ReadRouteConfig(layout);
		parser.ReadNet(layout);

		RoutingDB routingDB(designName, layout);

		routingDB.initGlobalNets(layout);

		routingDB.initEdges(layout);
		routingDB.initGlobalEdgeProfile(layout);

// Read routing file
		if (rtFile) {
			routingDB.readGlobalWires(layout, rtFile);
			routingDB.initRoutingTreeForAllGnets();
			routingDB.updateEdgeDemands();
			auto vpins = routingDB.getVpins(layout, layer - 1, twoCutNetsOnly, floatingNets, excludeNI);

			std::unordered_map<Gcell, size_t, HashGcell3d> uMapVpin;
			for (size_t i = 0; i < vpins.size(); i++) {
				auto &vp = vpins[i];
				uMapVpin[Gcell(vp.gnetID, vp.xCoord, vp.yCoord)] = i;
			}
            unordered_set<size_t> excludedGnetIds;
			for (size_t i = 0; i < vpins.size(); i++) {
				auto &vp = vpins[i];
				excludedGnetIds.insert(vp.gnetID);
			}
			if (outputRtFile) {
				cout << "Writing rtFile" << endl;
                size_t netCount = 0;
                if (!brokenNetsOnly) 
				    netCount = routingDB.writeGlobalWires(layout, outputRtFile, excludedGnetIds);
                for (size_t i = 0; i < vpins.size(); i++) {
                    auto &vp = vpins[i];
                    routingDB.writeGlobalWiresVpin(layout, outputRtFile, vp, netCount);
                }
			}
            if (outputNetsFile) {
            	routingDB.writeNets(layout, outputNetsFile, vpins, brokenNetsOnly);
            }
            if (outputKeyFile) {
				routingDB.writeKey(layout, outputKeyFile, vpins);            
            }
            
            if (obfuscate) {
                cout << "Start obfuscation, which may take a long time ..." << endl;
                int pass = 0, fail = 0;
                mt19937 gen;
                size_t maxTrials = 10;
                double stdevX = 0.01 * layout.getSizeX();
                double stdevY = 0.01 * layout.getSizeY();
                vector<size_t> randperm(vpins.size());
                iota(randperm.begin(), randperm.end(), 0);
                auto rng = mt19937{};
                shuffle(randperm.begin(), randperm.end(), rng);
                int i = 0;
                for (auto idx : randperm) {
                    cout << i << " / " << vpins.size() << endl;
                    i++;
                    auto res = routingDB.randomTraverse(layout, layer, vpins, idx, maxTrials, stdevX, stdevY, gen);
                    if (res.second) {
                        bool ret = routingDB.applyChange(layout, layer, vpins, idx, res.first);
                        assert(ret);
                        pass++;
                    } else {
                        fail++;
                    }
                }
                cout << "Vpin changed: " << pass << ", Not changed: " << fail << endl;
            }
			if (outputRtFile) {
                char obfusRtFile[100];
                sscanf(obfusRtFile, "%s.obfuscated.rt", outputRtFile);
				cout << "Writing obfuscated rtFile" << endl;
                size_t netCount = 0;
                if (!brokenNetsOnly) 
				    netCount = routingDB.writeGlobalWires(layout, obfusRtFile, excludedGnetIds);
                for (size_t i = 0; i < vpins.size(); i++) {
                    auto &vp = vpins[i];
                    routingDB.writeGlobalWiresVpin(layout, obfusRtFile, vp, netCount);
                }
			}
		} // if (rtFile)
    } // if (auxFile)
    return 0;
}
