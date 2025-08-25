#pragma once

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <map>
#include <tuple>

#include "Bump.hpp"
#include "DesignRule.hpp"

using namespace std;

void ParseBump(const string &inputPath, vector<Bump> &bumps, vector<double> &coordinates) ; // parse via_layer_X
void ParseDesignRule(const string &inputPath, DesignRule &designRule) ; // parse design rule
void ParseOffsetBump(const string &inputPath, vector<Bump> &bumps, vector<Bump> &matchedBumps) ; // parse offset_via_layer_X
void ParseNet(const string &inputPath, vector<Net> &nets) ; // parse netlist_X
void ParseTeardrop(const string &inputPath, vector<tuple<Bump, double, double, double, double>> &teardrops) ; // parse teardrop_X

// void parseOffestVia(vector<Bump>& offset_vias, string layer_to_offset_viaFile){
//     ifstream file(layer_to_offset_viaFile);
//     if (!file.is_open()) 
//         return;
    

//     string line;
//     while (getline(file, line)) {
//         istringstream iss(line);
//         string dieName, type;
//         int id;
//         double x, y, x2, y2;
//         if (!(iss >> dieName >> type >> id >> x >> y >> x2 >> y2)) {
//             cerr << "Offset via parse error：" << line << endl;
//             continue;
//         }
        
//         offset_vias.push_back(Bump(dieName, type, id, x, y));
//         offset_vias.push_back(Bump(dieName, type, id, x2, y2));
//     }
//     file.close();
// }

// // read design rule


// void parseTriangleEdge(vector<pair<double, double>>& triangleEdgeSource, vector<pair<double, double>>& triangleEdgeTarget, const string &fileName){
//     ifstream file(fileName);
//     string line;
//     while (getline(file, line)) {
//         istringstream iss(line);
//         double x1, y1, x2, y2;
//         // 解析每行的四個座標
//         if (iss >> x1 >> y1 >> x2 >> y2) {
//             triangleEdgeSource.push_back({x1, y1});
//             triangleEdgeTarget.push_back({x2, y2});
//         }
//     }

//     file.close();
// }
