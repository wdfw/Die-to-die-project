#include <iostream>
#include <vector>
#include <ctime>

#include "Bump.hpp"
#include "DesignRule.hpp"
#include "Parser.hpp"
#include "GlobalRouting.hpp"

//Drawer.hppㄧ定要include在最上層

using namespace std ; 

int main(int argc, char *argv[]) {

    // return 0 ;
    if (argc!=4) {
        cerr << "Usage: " << argv[0] << " <input_file> <design_rule_file> <output_dir>" << "\n";
        return 1 ;
    }

    string bumpFilePath = argv[1] ; 
    string designRulePath = argv[2] ;
    string outputDirectories = argv[3] ;

    DesignRule designRule ; 
    vector<Bump> allBumps, die1Bumps, die2Bumps ; 
    vector<double> coordinate ; 
    vector<RoutingGraph> allRDL ;
    vector<clock_t> globalRouteTimes, getailedRouteTimes ;

    ParseBump(bumpFilePath, allBumps, coordinate) ;
    ParseDesignRule(designRulePath, designRule) ;

    GlobalRoute(allBumps, coordinate, designRule, outputDirectories, allRDL, 0, globalRouteTimes);

    return 0 ;
}
// 目前假設情境
// 1. DIE1在DIE2左邊
// 2. DIE1與DIE2之的間距大於相鄰 bumps 最小X距離

//./bin/D2D case/d2d_case_bump.location case/design.rule d2d_result 48
//./bin/D2D case/d2d_case_bump.location case/design.rule wu_result/ 

//./bin/ShowResult result/ case/design.rule

//./bin/ShowResult wu_result/ case/design.rule
