#pragma once

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits> 
#include <iomanip> 
#include <filesystem>

#include "Bump.hpp"
#include "DesignRule.hpp"

using namespace std ; 

// generate Die‑1/Die‑2 bump coordinate files
// d2d_case_bump.location => 傳遞資訊到 bump layer 的點結構 => vector<Bump> bumps
// RDL1/via_layer_1 => 傳遞資訊到 RDL layer 的點結構 => vector<Bump> vias
void GenerateBumpCaseFiles(int numSignal, const string &casePath, const string &viaPath, double &offset) ;