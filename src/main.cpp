#include <iostream>

#include "Drawer.hpp"
#include "Bump.hpp"
#include "DesignRule.hpp"
#include "Parser.hpp"

//Drawer.hppㄧ定要include在最上層

using namespace std ; 
namespace bg = boost::geometry ;

int main(int argc, char *argv[]) {

    // return 0 ;
    if (argc!=4) {
        cerr << "Usage: " << argv[0] << " <input_file> <design_rule_file> <output_dir>" << "\n";
        return 1 ;
    }

    string bumpFilePath = argv[1] ; 
    string designRulePath = argv[2] ;
    string outputDirectories = argv[2] ;

    DesignRule designRule ; 
    vector<Bump> bumps ; 
    vector<double> coordinate ; 
    vector<RoutingGraph> allRDL ;
    vector<clock_t> globalRouteTimes, getailedRouteTimes ;

    ParseBump(bumpFilePath, bumps, coordinate) ;
    ParseDesignRule(designRulePath, designRule) ;
    
    num_of_layers = routingGuideGenerate(bumps, coordinate, designRule, outputDirectories, allRDL, offset, GlobalRouteTime);

    return 0 ;
}
//./bin/D2D case/d2d_case_bump.location case/design.rule d2d_result 48

//這裡的via N代表第N與N-1層的via