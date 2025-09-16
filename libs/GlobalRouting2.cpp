#include "GlobalRouting2.hpp"


//------------------------------------------ GlobalRouter Method Begin ------------------------------------------ 

GlobalRouter::GlobalRouter(const DesignRule& designRule, const vector<Bump>& bumps, const vector<double>& coordinate){
    this->designRule = designRule ; 
    this->bumps = bumps ; 
    this->coordinate = coordinate ; 
}

void GlobalRouter::Initial(){
    this->routingTimes.clear() ; 
    this->routingGraphs.clear() ; 
}

void GlobalRouter::GlobalRoute(const string& outputDirectories){
    double horizontalSpace, minY = numeric_limits<double>::max() ; 
    string directoryPath ; 

    Timer timer;

    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> marginalBumps ;
    
    //找出最左邊的Signal bumps
    for(const auto& bump : bumps){
        if(bump.type == SIGNAL){
            if(bump.y<minY){
                minY = bump.y ;
                marginalBumps.clear() ; 
            }
            if(fabs(bump.y - minY) < EPSILON_Y){
               marginalBumps.push_back(bump) ;
            }
        }
    }
    sort(marginalBumps.begin(), marginalBumps.end(), [](const Bump& a, const Bump& b) {return a.x < b.x;}) ;
    horizontalSpace = marginalBumps[1].x - marginalBumps[0].x ; 

    for(int layer = 1, leftBumpCount = bumps.size(); leftBumpCount; ++layer){
        directoryPath = outputDirectories + "RDL" + to_string(layer) + "/";

        if (!filesystem::exists(directoryPath)) filesystem::create_directories(directoryPath);
        
        timer.SetClock() ;

        ConstructRoutingGraph(bumps, routingCoordinate, designRule, directoryPath, horizontalSpace, layer, allRDL);
        globalRouteTimes.push_back(timer.GetDurationMilliseconds()); 

        leftBumpCount = 0 ;
        // for(const auto& bump : die1Bumps) if(bump.type==SIGNAL) ++leftBumpCount ;
    }
}


void GlobalRouter::FindLeftmostInEachRow(const vector<Bump>& bumps, vector<Bump>& leftMostBumps) {
    vector<vector<Bump>> groups ;
    vector<Bump> sorted_vias = bumps; sort(sorted_vias.begin(), sorted_vias.end(), [](const Bump& a, const Bump& b) {return a.y < b.y;});

    leftMostBumps.clear() ; 

    for (const auto& via : sorted_vias) {
        bool added = false;
        for (auto& group : groups) {
            if (!group.empty() && fabs(group[0].y - via.y) < epsilonY) {
                group.push_back(via);
                added = true;
                break;
            }
        }
        if (!added)  groups.push_back({via}) ;
    }

    // 每組選出 x 最小的 bump（即最左點）
    for(const auto& group : groups) leftMostBumps.push_back( *min_element(group.begin(), group.end(), [](const Bump& a, const Bump& b) {return a.x < b.x;}) );
}
void GlobalRouter::PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, double minimumHorizontalSpace){
    double min_x = coordinate[0], max_x = coordinate[2], current_x ;
    Bump dummyBump ;

    vector<Bump> leftmostBumps ; FindLeftmostInEachRow(bumps, leftmostBumps) ;
    
    paddingBumps = bumps ; 
    for (const auto& point : leftmostBumps) {
        double current_x = point.x;

        // 向左移動
        while (current_x - minimumHorizontalSpace >= min_x + minimumHorizontalSpace) {
            current_x -= minimumHorizontalSpace ;
            dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&dummyBump](const Bump& bump){return fabs(bump.x - dummyBump.x) < epsilonX && fabs(bump.y - dummyBump.y) < epsilonY;})==bumps.end()){
                paddingBumps.push_back(dummyBump);
            }
        }

        current_x = min_x ; 
        dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&dummyBump](const Bump& bump){return fabs(bump.x - dummyBump.x) < epsilonX && fabs(bump.y - dummyBump.y) < epsilonY;})==bumps.end()){
            paddingBumps.push_back(dummyBump);
        }

        // 向右移動
        current_x = point.x;
        while (current_x + minimumHorizontalSpace <= max_x - minimumHorizontalSpace) {
            current_x += minimumHorizontalSpace;
            dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&dummyBump](const Bump& bump){return fabs(bump.x - dummyBump.x) < epsilonX && fabs(bump.y - dummyBump.y) < epsilonY;})==bumps.end()){
                paddingBumps.push_back(dummyBump);
            }
        }

        current_x = max_x ; 
        dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&dummyBump](const Bump& bump){return fabs(bump.x - dummyBump.x) < epsilonX && fabs(bump.y - dummyBump.y) < epsilonY;})==bumps.end()){
            paddingBumps.push_back(dummyBump);
        }
    }
}

void GlobalRouter::ConstructRoutingGraph(const string& outputPath, int layer, double minimumHorizontalSpace){
    RoutingGraph RDL1, RDL2 ; //RDL1 for left part, RDL2 for right part
    
    vector<Bump> die1Bumps, die2Bumps, die1PaddingBumps, die2PaddingBumps ; 
    vector<double> die1Coordinate, die2Coordinate ; 

    ofstream outputViaFile, outputTriangulationFile, outputNetlistFile ; 
    die1Coordinate = routingCoordinate ; coordinate[2] = coordinate[2]/2 - minimumHorizontalSpace ; 
    die2Coordinate = routingCoordinate ; coordinate[0] = coordinate[2]/2 + Hor_SPACINminimumHorizontalSpaceG_X ; 

    copy_if(bumps.begin(), bumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(bumps.begin(), bumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;


    PaddingBumps(die1Bumps, die1PaddingBumps, minimumHorizontalSpace) ;
    PaddingBumps(die2Bumps, die2PaddingBumps, minimumHorizontalSpace) ;

    for(auto bump : die1Bumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    for(auto bump : die1PaddingBumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    for(auto bump : die2Bumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    for(auto bump : die2PaddingBumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    // // 使用 CDT 進行三角化
    // RDL1.edge_nodes = Triangulation(RDL1.via_nodes, designRule);
    // RDL2.edge_nodes = Triangulation(RDL2.via_nodes, designRule);
    
    // FindAdjacentViaNodes(RDL1) ; //!!!!!
    // FindAdjacentViaNodes(RDL2) ; //!!!!!

    // AddAccessViaEdges(RDL1) ;
    // AddAccessViaEdges(RDL2) ;

    // AddCrossTileEdges(RDL1, designRule) ;
    // AddCrossTileEdges2(RDL1, designRule) ;

    outputViaFile.open(output_path + "via_layer_" + to_string(layer)) ;
    for(auto& bump : die1Bumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die2Bumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die1DummyBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die2DummyBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : debugBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    outputViaFile.close() ;

    outputTriangulationFile.open(output_path + "triangulation_edge");
    for(auto& edgeNode : RDL1.edge_nodes) outputTriangulationFile << edgeNode.start.first << " " << edgeNode.start.second << " " << edgeNode.end.first << " " << edgeNode.end.second << "\n" ;
    for(auto& edgeNode : RDL2.edge_nodes) outputTriangulationFile << edgeNode.start.first << " " << edgeNode.start.second << " " << edgeNode.end.first << " " << edgeNode.end.second << "\n" ;
    for(auto& edgeNode : debugEdges) outputTriangulationFile << edgeNode.start.first << " " << edgeNode.start.second << " " << edgeNode.end.first << " " << edgeNode.end.second << "\n" ;
    outputTriangulationFile.close() ;

    outputNetlistFile.open(output_path + "netlist_" + to_string(layer));
    for(auto& net : debugNets){
        for(auto& seg : net){
            outputNetlistFile << net.name << " " << get<0>(seg) << " " << get<1>(seg) << " " << get<2>(seg) << " " << get<3>(seg) << "\n" ;
        }
    }
    outputNetlistFile.close() ;
}

//------------------------------------------ GlobalRouter Method End ------------------------------------------ 
