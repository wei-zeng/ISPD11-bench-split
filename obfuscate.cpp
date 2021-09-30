/*
 * obfuscate.cpp
 *
 *  Created on: Aug 28, 2019
 *      Author: wzeng
 */
#include "parser.h"
#include "grDB.h"
#include <chrono>
#include <set>
#include <cassert>

std::string ObfusOption::toString() const {
	std::stringstream ss;
	ss << "[" << gnetID << "," << d0 << "," << d1 << "]";
	return ss.str();
}

std::pair<ObfusOption, bool> RoutingDB::randomTraverse(const Layout &layout,
		int layer, vector<Vpin> &vps, size_t best_vpin, size_t maxTrials,
		double stdevX, double stdevY, std::mt19937 &rng) {
	vector<ObfusOption> opts;
	RoutingDB &tmpRtDB = (*this);
	int gnetID = vps[best_vpin].gnetID;
	Gnet cleanGnet = _gnets[gnetID];
	assert(cleanGnet.isGlobal());
	Vpin vp = vps[best_vpin];
	auto res = tmpRtDB.ripUpPin(layout, vp, layer - 1);
	Gnet rippedUpGnet = tmpRtDB._gnets[gnetID];
	std::normal_distribution<double> randnX(0.0, stdevX);
	std::normal_distribution<double> randnY(0.0, stdevY);

    int trial;
	for (trial = 0; trial < maxTrials; trial++) {
		auto vp0 = Vpin(vp);
		double rnd0 = randnX(rng);
        double rnd1 = randnY(rng);
        int d0 = round(rnd0 / layout._cellWidth);
		int d1 = round(rnd1 / layout._cellHeight);
        cout << "Trial " << trial << " d0 = " << d0 << " d1 = " << d1 << endl;
        int margin0, margin2;
		vp0.xCoord += d0;
        vp0.yCoord += d1;
        if (vp0.xCoord < 0 || vp0.xCoord > static_cast<int>(layout._numTilesX) - 1)
            continue;
		if (vp0.yCoord < 0 || vp0.yCoord > static_cast<int>(layout._numTilesY) - 1)
			continue;
        if (d0 != 0 || d1 != 0) {
            auto vp1 = Vpin(vps[vp0.matchingVpinIdx]);
            auto res2 = tmpRtDB.rerouteNet(layout, vp0,
                    vps[vp0.matchingVpinIdx].gnetGrids, 0, layer - 1, 15, false,
                    false, false);
            if (std::get<0>(res2) != INFINITY) {
                margin0 = std::get<3>(res2);
            } else {
                continue;
            }
            //	cout << "(Reroute) Failed!" << endl;
            auto res3 = tmpRtDB.rerouteNet(layout, vp0, vps[vp0.matchingVpinIdx],
                    layer, 8, 15, false, false);

            if (std::get<0>(res3) != INFINITY) {
                margin2 = std::get<3>(res3);
            } else {
                continue;
            }
        }

		ObfusOption opt;
		opt.gnetID = gnetID;
		opt.d0 = d0;
		opt.d1 = d1;
		opt.margin0 = margin0;
		opt.margin2 = margin2;
		opts.push_back(opt);
		if (tmpRtDB.isDirty()) {
			tmpRtDB.restoreNet(gnetID, rippedUpGnet, _gnets[gnetID]._ruEdges,
					_gnets[gnetID]._ruVias, true);
			tmpRtDB.setDirty(); // not restored to clean yet.
		}
		break;
	} // for trial
    if (trial == maxTrials) {
        cout << "No feasible move for this v-pin." << endl;
    }
	tmpRtDB.restoreNet(gnetID, cleanGnet, _gnets[gnetID]._cleanEdges,
		_gnets[gnetID]._cleanVias, true);
	assert(_gnets[gnetID] == cleanGnet);
	for (auto opt : opts) {
		cout << "Found feasible move: " << opt.toString() << endl;
	}
	return opts.size() > 0 ?
			make_pair(opts[0], true) : make_pair(ObfusOption(), false);
}

bool RoutingDB::applyChange(const Layout &layout,
		int layer, vector<Vpin> &vps, size_t best_vpin, const ObfusOption &opt) {
    if (opt.d0 == 0 && opt.d1 == 0) {
        return true;
    }
	int gnetID = vps[best_vpin].gnetID;
    assert(gnetID == opt.gnetID);
	Gnet cleanGnet = _gnets[gnetID];
	assert(cleanGnet.isGlobal());
	Vpin &vp0 = vps[best_vpin];
	auto res = ripUpPin(layout, vp0, layer - 1);
	int margin0 = opt.margin0;
	int margin2 = opt.margin2;
	vp0.xCoord += opt.d0;
    vp0.yCoord += opt.d1;
    if (vp0.xCoord < 0 || vp0.xCoord > static_cast<int>(layout._numTilesX) - 1) {
        restoreNet(gnetID, cleanGnet, _gnets[gnetID]._cleanEdges, 
                    _gnets[gnetID]._cleanVias, true);
        return false;
    }
    if (vp0.yCoord < 0 || vp0.yCoord > static_cast<int>(layout._numTilesY) - 1) {
        restoreNet(gnetID, cleanGnet, _gnets[gnetID]._cleanEdges, 
                    _gnets[gnetID]._cleanVias, true);
        return false;
    }
    auto vp1 = Vpin(vps[vp0.matchingVpinIdx]);
    auto res2 = rerouteNet(layout, vp0,
            vp1.gnetGrids, 0, layer - 1, opt.margin0, true,
            true, !vp0.hasOriginalRoot);
    if (std::get<0>(res2) == INFINITY) {
        restoreNet(gnetID, cleanGnet, _gnets[gnetID]._cleanEdges, 
                    _gnets[gnetID]._cleanVias, true);
        return false;
    }
    //	cout << "(Reroute) Failed!" << endl;
    auto res3 = rerouteNet(layout, vp0, vp1,
            layer, 8, opt.margin2, true, true);

    if (std::get<0>(res3) == INFINITY) {
        restoreNet(gnetID, cleanGnet, _gnets[gnetID]._cleanEdges, 
                    _gnets[gnetID]._cleanVias, true);
        return false;
    }
    return true;
}
