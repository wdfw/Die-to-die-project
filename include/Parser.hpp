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
#include "Utils.hpp"

using namespace std;

void ParseBump(const string &inputPath, vector<Bump> &bumps, vector<double> &coordinates) ; // 爬Bump檔, 檔名為via_layer_X
void ParseDesignRule(const string &inputPath, DesignRule &designRule) ; // 爬Design rule檔
void ParseOffsetBump(const string &inputPath, vector<Bump> &bumps, vector<Bump> &matchedBumps) ; // 爬Offset bump檔(包含bump與其offset via的位置), 檔名為offset_via_layer_X
void ParseNet(const string &inputPath, vector<Net> &nets) ; // 爬Net檔, 檔名為netlist_X
void ParseTeardrop(const string &inputPath, vector<tuple<Bump, double, double, double, double>> &teardrops) ; // 爬Teardrop檔, 檔名為tear_drop_X, (註:因為teardrop在演算法中用不到, 所以直接用tuple來表示)
void ParseTriangulation(const string &inputPath, vector<tuple<double, double, double, double>> &triangulationEdges) ; // 爬Triangulation檔, 檔名為triangulation_edge
 
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
