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
    vector<Bump> allBumps ; 
    vector<double> routingCoordinate ; 
    vector<RoutingGraph> allRDL ;
    vector<clock_t> globalRouteTimes, getailedRouteTimes ;

    ParseBump(bumpFilePath, allBumps, routingCoordinate) ;
    ParseDesignRule(designRulePath, designRule) ;

    GlobalRoute(allBumps, routingCoordinate, designRule, outputDirectories, allRDL, 0, globalRouteTimes);

    return 0 ;
}
// 目前假設情境
// 1. DIE1在DIE2左邊
// 2. DIE1與DIE2之的間距大於相鄰 bumps 最小X距離

//./bin/D2D case/d2d_case_bump.location case/design.rule d2d_result 48 12 #YY的版本 D2D  bump_file design_rule_file result_folder/ #bump offset
//./bin/D2D case/d2d_case_bump.location case/design.rule wu_result/  #我的版本 D2D bump_file design_rule_file result_folder/

//./bin/ShowResult result/ case/design.rule #根據design_rule顯示folder的資料
//./bin/ShowResult wu_result/ case/design.rule

// ./bin/D2D result/RDL1/via_layer_1 case/design.rule wu_result/ > wu_result/RDL1/via_layer_1 ; ./bin/ShowResult wu_result/ case/design.rule