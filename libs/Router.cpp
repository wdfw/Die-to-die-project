#include "Router.hpp"

vector<Bump> debugBumps ; 
// vector<EdgeNode> debugEdges ; 
vector<Net> debugNets ; 
vector<tuple<string, double, double>> debugLabels ; 
//------------------------------------------ Router Method Begin ------------------------------------------ 

Router::Router(const DesignRule& designRule, const vector<Bump>& bumps, const vector<double>& coordinate){
    this->designRule = designRule ; 
    this->bumps = bumps ; 
    this->coordinate = coordinate ; 
}

void Router::Initialize(){
    this->routingTimes.clear() ; 
    this->routingGraphs.clear() ; 
}
void Router::FindHorizontalSpace(const vector<Bump>& bumps, double& horizontalSpace){
    vector<vector<Bump>> groups ;
    horizontalSpace = numeric_limits<double>::max() ;

    Matrixization(bumps, groups, epsilonY) ; 
    for(int i=0; i<groups.size(); ++i){
        for(int j=0; j<groups[i].size()-1; ++j){
            horizontalSpace = min(horizontalSpace, groups[i][j+1].x - groups[i][j].x) ; 
        }
    }
}
void Router::SelectRoutingBumps(const vector<Bump>& bumps, vector<Bump>& routingBumps, vector<Bump>& offsetBumps){
    int leftSignalCount = 0 ;
    set<int> selectedID = {36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47} ; 
    routingBumps.clear() ; 
    offsetBumps.clear() ; 
    for(int i=0; i<bumps.size(); i++){
        if(bumps[i].type==VSS){
            routingBumps.push_back(bumps[i]) ;
            offsetBumps.push_back(bumps[i]) ;
        }else if(bumps[i].type==VDD){
            routingBumps.push_back(bumps[i]) ;
        }else if(bumps[i].type==SIGNAL){
            if(selectedID.find(bumps[i].id)!=selectedID.end()){
                routingBumps.push_back(bumps[i]) ;
            }else{
                offsetBumps.push_back(bumps[i]) ;
                ++ leftSignalCount ; 
            }
        }
    }
    
    if(!leftSignalCount) offsetBumps.clear() ; 
}

void Router::Solve(const string& outputDirectories){
    double horizontalSpace, minY = numeric_limits<double>::max() ; 
    string directoryPath ; 

    Timer timer;

    vector<Bump> unroutedBumps = bumps ;
    vector<Bump> marginalBumps, routingBumps, offsetBumps, viaBumps ;
    
    for(int layer = 1, leftBumpCount = bumps.size(); leftBumpCount; ++layer){
        RoutingGraph2 routingGraph ; 
        directoryPath = outputDirectories + "RDL" + to_string(layer) + "/";

        if (!filesystem::exists(directoryPath)) filesystem::create_directories(directoryPath);
        
        timer.SetClock() ;

        FindHorizontalSpace(unroutedBumps, horizontalSpace) ;
        SelectRoutingBumps(unroutedBumps, routingBumps, offsetBumps) ;
        CreateViaBumps(offsetBumps, viaBumps) ;
        ConstructRoutingGraph(routingBumps, offsetBumps, viaBumps, routingGraph, horizontalSpace) ;

        GlobalRoute(routingBumps, routingGraph) ; 
        
        GenerateGraphFile(routingBumps, offsetBumps, viaBumps, routingGraph, layer, directoryPath) ;

        unroutedBumps = offsetBumps ; 
        routingTimes.push_back(timer.GetDurationMilliseconds()); 

        leftBumpCount = 0 ; // for(const auto& bump : die1Bumps) if(bump.type==SIGNAL) ++leftBumpCount ;
    }
}

void Router::CreateViaBumps(const vector<Bump>& offsetBumps, vector<Bump>& viaBumps){
    for(auto& bump : offsetBumps){
        viaBumps.push_back(bump) ; 
        viaBumps.back().x = bump.x -( designRule.minimumViaSpacing + designRule.viaOpeningDiameter) ;
        viaBumps.back().y = bump.y ;
    }
}

void Router::GenerateGraphFile(const vector<Bump>& routingBumps, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps, const RoutingGraph2& graph, int layer,  const string& directoryPath){
    ofstream outputViaFile, outputOffsetViaFile, outputTriangulationFile, outputNetlistFile, outputDebugLabelFile ; 

    outputViaFile.open(directoryPath + "via_layer_" + to_string(layer)) ;
    for(auto& viaNode : graph.viaNodes) outputViaFile << Bump2Str(*viaNode) << "\n" ;
    // for(auto& bump : debugBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    outputViaFile.close() ;

    outputOffsetViaFile.open(directoryPath + "offset_via_layer_" + to_string(layer)) ;
    for(int i=0; i<offsetBumps.size(); i++){
        outputOffsetViaFile << Bump2Str(offsetBumps[i]) << " " << viaBumps[i].x << " " << viaBumps[i].y << "\n" ;
    }
    outputOffsetViaFile.close() ;

    outputTriangulationFile.open(directoryPath + "triangulation_edge");
    for(auto& RDL : {graph}){
        for(auto& start : RDL.viaNodes){
            for(auto& traget : start->viaNodes){
                outputTriangulationFile << start->x << " " << start->y << " " << traget->x << " " << traget->y << "\n" ;
            }
        }
        // for(auto& start : RDL.tileNodes){
        //     for(auto& traget : start->tileNodes){
        //         outputTriangulationFile << traget.crossedViaNode1->x << " " << traget.crossedViaNode1->y << " " << traget.crossedViaNode2->x << " " << traget.crossedViaNode2->y << "\n" ;
        //     }
        // }
        for(auto& start : RDL.tileNodes){
            for(auto& traget : start->tileNodes){
                outputTriangulationFile << start->x << " " << start->y << " " << traget->x << " " << traget->y << "\n" ;
            }
        }
    }
    // for(auto& viaNode : debugEdges) outputTriangulationFile << edgeNode.start.first << " " << edgeNode.start.second << " " << edgeNode.end.first << " " << edgeNode.end.second << "\n" ;
    outputTriangulationFile.close() ;
//equal length

    outputNetlistFile.open(directoryPath + "netlist_" + to_string(layer));
    for(auto& net : debugNets){
        for(auto& seg : net){
            outputNetlistFile << net.name << " " << get<0>(seg) << " " << get<1>(seg) << " " << get<2>(seg) << " " << get<3>(seg) << "\n" ;
        }
    }
    outputNetlistFile.close() ;

    outputDebugLabelFile.open(directoryPath + "debug_label");
    for(auto& debugLabel : debugLabels){
        outputDebugLabelFile <<  get<0>(debugLabel) << " " << get<1>(debugLabel) << " " << get<2>(debugLabel) << "\n" ;
    }
    outputDebugLabelFile.close() ;
}

void Router::PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, const vector<double>& coordinate, double minimumHorizontalSpace){
    double min_x = coordinate[0], max_x = coordinate[2], current_x ;
    Bump dummyBump ;

    vector<Bump> leftmostBumps ; FindLeftMostInEachRow(bumps, leftmostBumps) ;
    // for(auto bump:leftmostBumps) debugBumps.push_back(bump) ;
    paddingBumps.clear() ;
    for (const auto& point : leftmostBumps) {
        double current_x = point.x;

        // 向左移動
        while (current_x - minimumHorizontalSpace >= min_x + minimumHorizontalSpace) {
            current_x -= minimumHorizontalSpace ;
            dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);
            
            if(find_if(bumps.begin(), bumps.end(), [&dummyBump, this](const Bump& bump){return fabs(bump.x - dummyBump.x) < this->epsilonX && fabs(bump.y - dummyBump.y) < this->epsilonY;})==bumps.end()){
                paddingBumps.push_back(dummyBump);
            }
        }

        current_x = min_x ; 
        dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&dummyBump,this](const Bump& bump){return fabs(bump.x - dummyBump.x) < this->epsilonX && fabs(bump.y - dummyBump.y) < this->epsilonY;})==bumps.end()){
            paddingBumps.push_back(dummyBump);
        }

        // 向右移動
        current_x = point.x;
        while (current_x + minimumHorizontalSpace <= max_x - minimumHorizontalSpace) {
            current_x += minimumHorizontalSpace;
            dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&dummyBump, this](const Bump& bump){return fabs(bump.x - dummyBump.x) < this->epsilonX && fabs(bump.y - dummyBump.y) < this->epsilonY;})==bumps.end()){
                paddingBumps.push_back(dummyBump);
            }
        }

        current_x = max_x ; 
        dummyBump = Bump("Dummy", DUMMY, paddingBumps.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&dummyBump, this](const Bump& bump){return fabs(bump.x - dummyBump.x) < this->epsilonX && fabs(bump.y - dummyBump.y) < this->epsilonY;})==bumps.end()){
            paddingBumps.push_back(dummyBump);
        }
       
    }
}


void Router::Triangulation(RoutingGraph2& graph) {
    CDT cdt;
   
    vector<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>> triangularEdges  ;
    vector<Point> allPoints;
    vector<vector<Point>> y_groups; // === Y 分組（允許誤差）===
    map<Point, shared_ptr<ViaNode2>> point_to_via_ptr;

    for (const auto &node : graph.viaNodes) {
        Point p(node->x, node->y);
        point_to_via_ptr[p] = node;
        allPoints.push_back(p);
    }

    sort(allPoints.begin(), allPoints.end(), [](const Point& a, const Point& b) {return a.y() < b.y();});

    for (const auto& pt : allPoints) {
        bool added = false;
        for (auto& group : y_groups) {
            if (!group.empty() && fabs(group[0].y() - pt.y()) < epsilonY) {
                group.push_back(pt);
                added = true;
                break;
            }
        }
        if (!added) y_groups.push_back({pt});
    }

    // === 每組內按 x 連線 ===
    for (const auto& group : y_groups) {
        vector<Point> sorted_group = group;
        sort(sorted_group.begin(), sorted_group.end(), [](const Point& a, const Point& b) {
            return a.x() < b.x();
        });
        for (size_t i = 1; i < sorted_group.size(); ++i) { // 我的 CDT 規則 => 相同 y 值的點, 要先連在一起
            cdt.insert_constraint(sorted_group[i - 1], sorted_group[i]);
        }
    }

    // === 每組最左點/最右點垂直連線 === => 我的 CDT 規則
    vector<Point> left_points, right_points;
    for (const auto& group : y_groups) {
        auto [min_it, max_it] = minmax_element(group.begin(), group.end(),
            [](const Point& a, const Point& b) { return a.x() < b.x(); });
        left_points.push_back(*min_it);
        right_points.push_back(*max_it);
    }

    sort(left_points.begin(), left_points.end(), [](const Point& a, const Point& b) {
        return a.y() < b.y();
    });
    sort(right_points.begin(), right_points.end(), [](const Point& a, const Point& b) {
        return a.y() < b.y();
    });

    for (size_t i = 1; i < left_points.size(); ++i)
        cdt.insert_constraint(left_points[i - 1], left_points[i]);
    for (size_t i = 1; i < right_points.size(); ++i)
        cdt.insert_constraint(right_points[i - 1], right_points[i]);

    // === 插入所有點三角化 ===
    cdt.insert(allPoints.begin(), allPoints.end());

    // === 處理 CDT 邊 ===
    int edgeId = 0;
    
    // ofstream outFile(output_path + "triangulation_edge", ofstream::app);

    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        auto segment = cdt.segment(*edge);
        double dx = segment.source().x() - segment.target().x();
        double dy = segment.source().y() - segment.target().y();
        double L = sqrt(dx * dx + dy * dy);
        double ux = dx / L, uy = dy / L;

        double newSourceX, newSourceY, newTargetX, newTargetY;
        if (point_to_via_ptr[segment.source()]->type == DUMMY) {
            newSourceX = segment.source().x();
            newSourceY = segment.source().y();
            // outFile << newSourceX << " " << newSourceY << " ";
        } else { // 因為CDT的邊, 都是從 bump 的圓心出發, 畫出來的邊會壓到圓形, 所以要特殊處理
            newSourceX = segment.source().x() - (1 + designRule.viaOpeningDiameter / 2 +
                         (designRule.viaPadDiameter - designRule.viaOpeningDiameter) / 2) * ux;
            newSourceY = segment.source().y() - (1 + designRule.viaOpeningDiameter / 2 +
                         (designRule.viaPadDiameter - designRule.viaOpeningDiameter) / 2) * uy;
            // outFile << newSourceX << " " << newSourceY << " ";
        }

        if (point_to_via_ptr[segment.target()]->type == DUMMY) {
            newTargetX = segment.target().x();
            newTargetY = segment.target().y();
            // outFile << newTargetX << " " << newTargetY << "\n";
        } else { // 因為CDT的邊, 都是從 bump 的圓心出發, 畫出來的邊會壓到圓形, 所以要特殊處理
            newTargetX = segment.target().x() + (1 + designRule.viaOpeningDiameter / 2 +
                         (designRule.viaPadDiameter - designRule.viaOpeningDiameter) / 2) * ux;
            newTargetY = segment.target().y() + (1 + designRule.viaOpeningDiameter / 2 +
                         (designRule.viaPadDiameter - designRule.viaOpeningDiameter) / 2) * uy;
            // outFile << newTargetX << " " << newTargetY << "\n";
        }

        double x = (segment.source().x() + segment.target().x()) / 2; // CDT線段的中心點, 當作 EdgeNode 的 x, y 座標
        double y = (segment.source().y() + segment.target().y()) / 2;
        double distance = hypot(newSourceX - newTargetX, newSourceY - newTargetY);
        int capacity = floor(distance / (2*(designRule.minimumLineWidth + designRule.minimumLineSpacing))); //!!!!!
        // if(capacity > 5 || capacity < 5) capacity = 5; // 限制最大容量為 5

        point_to_via_ptr[segment.source()]->viaNodes.push_back(point_to_via_ptr[segment.target()]) ;
        point_to_via_ptr[segment.target()]->viaNodes.push_back(point_to_via_ptr[segment.source()]) ;

        triangularEdges.push_back({point_to_via_ptr[segment.source()], point_to_via_ptr[segment.target()]});
    }

    int n = triangularEdges.size() ; 

    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            for(int k=j+1; k<n; k++){
                map<shared_ptr<ViaNode2>, int> countPtr ; 
                double aveX = 0.0, aveY = 0.0 ;
                bool isClosed = true ;
                
                for(auto w : {i,j,k}){
                    ++ countPtr[triangularEdges[w].first] ;
                    ++ countPtr[triangularEdges[w].second] ;

                    aveX += triangularEdges[w].first->x ; aveX += triangularEdges[w].second->x ; 
                    aveY += triangularEdges[w].first->y ; aveY += triangularEdges[w].second->y ; 
                }
                

                for(auto& [v, c] : countPtr){
                    if(c!=2) isClosed = false ; 
                }
              
                if(isClosed){   
                    // debugBumps.emplace_back("DUMMY", DUMMY, 0, aveX/6, aveY/6) ;
                    // vector<ViaNode2> vias ; 
                    // for(auto& [c, v] : uniqueVias) vias.push_back(v) ;

                    // vector<Edge> edges = {
                    //     {RDL.edge_nodes[i].id, RDL.edge_nodes[i].vias, RDL.edge_nodes[i].capacity},
                    //     {RDL.edge_nodes[j].id, RDL.edge_nodes[j].vias, RDL.edge_nodes[j].capacity},
                    //     {RDL.edge_nodes[k].id, RDL.edge_nodes[k].vias, RDL.edge_nodes[k].capacity},
                    // } ;
                    shared_ptr<TileNode2> newTileNode = make_shared<TileNode2>("DUMMY", DUMMY, graph.tileNodes.size(), aveX/6, aveY/6) ; 

                    for(auto& [viaNode, c] : countPtr){
                        newTileNode->viaNodes.push_back(viaNode) ;
                        viaNode->tileNodes.push_back(newTileNode) ; 
                    }
                    
                    graph.tileNodes.push_back(newTileNode) ;
                }
            }
        }
    }
}

void Router::ConnectTileTileEdges(RoutingGraph2& graph){
    int n = graph.tileNodes.size() ; 

    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            map<shared_ptr<ViaNode2>, int> uniqueViaNodes ; 
            vector<shared_ptr<ViaNode2>> crossedViaNodes ; 
            for(auto& viaNode : graph.tileNodes[i]->viaNodes) ++ uniqueViaNodes[viaNode] ; 
            for(auto& viaNode : graph.tileNodes[j]->viaNodes) ++ uniqueViaNodes[viaNode] ; 

            if(uniqueViaNodes.size()==4){
                for(auto& [viaNode, count] : uniqueViaNodes){
                    if(count==2) crossedViaNodes.push_back(viaNode) ;
                }
                graph.tileNodes[i]->tileNodes.push_back(graph.tileNodes[j]) ;
                graph.tileNodes[j]->tileNodes.push_back(graph.tileNodes[i]) ;

                graph.tileNodes[i]->tileNodes.back().crossedViaNode1 = crossedViaNodes[0] ;
                graph.tileNodes[i]->tileNodes.back().crossedViaNode2 = crossedViaNodes[1] ;

                graph.tileNodes[j]->tileNodes.back().crossedViaNode1 = crossedViaNodes[0] ;
                graph.tileNodes[j]->tileNodes.back().crossedViaNode2 = crossedViaNodes[1] ;

                graph.tileNodes[j]->tileNodes.back().capacity = graph.tileNodes[i]->tileNodes.back().capacity ; 
            }
        }
    }
}

void Router::CreateViaNodes(RoutingGraph2& graph, const vector<Bump>& bumps){
    for(auto& bump : bumps){
        graph.viaNodes.push_back(make_shared<ViaNode2>(bump)) ;
    }
}
double cross(pair<double, double> o, pair<double, double> a, pair<double, double> b){
    return (a.first - o.first) * (b.second - o.second) - (a.second - o.second) * (b.first - o.first) ;
}

void Router::CombineRDLs(RoutingGraph2& graph1, RoutingGraph2& graph2, RoutingGraph2& graph3, const vector<double>& coordinate1, const vector<double>& coordinate2){
    double cx = (coordinate1[2] + coordinate2[0])/2 ;
    double cy = (coordinate1[1] + coordinate1[3] + coordinate2[1] + coordinate2[3])/4 ;

    vector<shared_ptr<TileNode2>> die1EscapeTileNodes, die2EscapeTileNodes ; 
    for(int i=0; i<2; i++){
        RoutingGraph2& graph = (i==0) ? graph1 : graph2 ; 
        vector<shared_ptr<TileNode2>>& dieEscapeTileNodes = (i==0) ? die1EscapeTileNodes : die2EscapeTileNodes ; 

        int n = graph.tileNodes.size(), m = graph.viaNodes.size() ; 

        for(int i=0; i<n; i++){
            int crossedCount = 0 ; 
            shared_ptr<TileNode2> rowTargetNode = make_shared<TileNode2>("Dummy", DUMMY, 0, cx, graph.tileNodes[i]->y) ; 
            vector<shared_ptr<ViaNode2>> crossedViaNodes ; 
            for(int j=0; j<m; j++){
                for(int k=0; k<graph.viaNodes[j]->viaNodes.size(); k++){
                    pair<double, double> pt1 = {rowTargetNode->x, rowTargetNode->y} ;
                    pair<double, double> pt2 = {graph.tileNodes[i]->x, graph.tileNodes[i]->y} ;
                    pair<double, double> pt3 = {graph.viaNodes[j]->x, graph.viaNodes[j]->y} ;
                    pair<double, double> pt4 = {graph.viaNodes[j]->viaNodes[k]->x, graph.viaNodes[j]->viaNodes[k]->y} ;
                    double d1 = cross(pt1, pt2, pt3) ;
                    double d2 = cross(pt1, pt2, pt4) ; 
                    double d3 = cross(pt3, pt4, pt1) ;
                    double d4 = cross(pt3, pt4, pt2) ;
                    if( ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) or (d3 < 0 && d4 > 0)) ){
                        ++crossedCount ;
                        crossedViaNodes = {graph.viaNodes[j], graph.viaNodes[j]->viaNodes[k]} ;
                    } 
                }
            }

            if(crossedCount<=2){
                dieEscapeTileNodes.push_back(graph.tileNodes[i]) ;
                dieEscapeTileNodes.back()->tileNodes.push_back(make_shared<TileNode2>("DUMMY", DUMMY, 0, 0, 0)) ;

                dieEscapeTileNodes.back()->tileNodes.back().crossedViaNode1 = crossedViaNodes[0] ;
                dieEscapeTileNodes.back()->tileNodes.back().crossedViaNode2 = crossedViaNodes[1] ;
            }
        }
    }

    sort(die1EscapeTileNodes.begin(), die1EscapeTileNodes.end(), [](const shared_ptr<TileNode2>& a, const shared_ptr<TileNode2>& b){return a->y<b->y; }) ;
    sort(die2EscapeTileNodes.begin(), die2EscapeTileNodes.end(), [](const shared_ptr<TileNode2>& a, const shared_ptr<TileNode2>& b){return a->y<b->y; }) ;

    if(die1EscapeTileNodes.size()!=die2EscapeTileNodes.size()) throw runtime_error("Invaild size on dies escape tile nodes (" + to_string(die1EscapeTileNodes.size()) + "," + to_string(die2EscapeTileNodes.size()) + ")") ; 

    // cout << die1EscapeTileNodes.size() << "\n" ;
    // cout << &die1EscapeTileNodes << " " << die1EscapeTileNodes.back()->tileNodes.back().crossedViaNode1 << "||\n" ;
    shared_ptr<ViaNode2> crossedViaNode1, crossedViaNode2 ; 

    for(int i=0; i<die1EscapeTileNodes.size(); i++){
        crossedViaNode1 = die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 ;
        crossedViaNode2 = die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 ;

        die1EscapeTileNodes[i]->tileNodes.back() = die2EscapeTileNodes[i] ;
        die2EscapeTileNodes[i]->tileNodes.back() = die1EscapeTileNodes[i] ;

        die1EscapeTileNodes[i]->tileNodes.back().capacity = die2EscapeTileNodes[i]->tileNodes.back().capacity ;

        die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 = crossedViaNode1 ; 
        die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 = crossedViaNode2 ; 

        die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 = crossedViaNode1 ; 
        die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 = crossedViaNode2 ; 

        // double  dx = crossedViaNode1->x - crossedViaNode2->x ; 
        // double dy = crossedViaNode1->y - crossedViaNode2->y ; 

        // cout <<  sqrt(dx*dx + dy*dy) << " " << designRule.viaPadDiameter << " " << designRule.minimumLineWidth << " " << designRule.minimumLineSpacing << "\n" ;
        // cout << crossedViaNode1->x << " " << crossedViaNode1->y << " " << crossedViaNode2->x << " " << crossedViaNode2->y << "\n\n" ; 
        // debugNets.push_back(Net("DUMMY", {{die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1->x, die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1->y, 
        //         die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2->x, die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2->y}})) ;
    }


    graph3 = RoutingGraph2(graph1) ; 
    for(auto& viaNode : graph2.viaNodes) graph3.viaNodes.push_back(viaNode) ;
    for(auto& tileNode : graph2.tileNodes) graph3.tileNodes.push_back(tileNode) ;

}

void Router::SetCapacity(RoutingGraph2& graph, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps){
    double dx, dy, ddx, ddy, distance, r1, r2 ;
    int overlapBumpIndex ; 

    for(auto& viaNode : graph.viaNodes){
        for(auto& tileNode : viaNode->tileNodes){
            (*tileNode.capacity) = numeric_limits<int>::max();
        }
    }
    for(auto& tileNode1 : graph.tileNodes){
        for(auto& tileNode2 : tileNode1->tileNodes){
            
            bool isOVerlap = false ;
            for(int i=0; i<viaBumps.size(); i++){
                const Bump& viaBump = viaBumps[i] ;
                bool cond1 = viaBump.type!=DUMMY ;
                bool cond2 = min(tileNode2.crossedViaNode1->x, tileNode2.crossedViaNode2->x) < viaBump.x && viaBump.x < max(tileNode2.crossedViaNode1->x, tileNode2.crossedViaNode2->x) ;
                bool cond3 = viaBump.y - designRule.viaPadDiameter/2 <= tileNode2.crossedViaNode1->y && tileNode2.crossedViaNode1->y <= viaBump.y + designRule.viaPadDiameter/2 ;
                bool cond4 = viaBump.y - designRule.viaPadDiameter/2 <= tileNode2.crossedViaNode2->y && tileNode2.crossedViaNode2->y <= viaBump.y + designRule.viaPadDiameter/2 ;

                if(cond1 && cond2 && cond3 && cond4){
                    isOVerlap = true ; overlapBumpIndex = i ; 
                    break ; 
                }
            }

            if(!tileNode2.crossedViaNode1 || !tileNode2.crossedViaNode2) continue ;

            if(isOVerlap){
                ddx = offsetBumps[overlapBumpIndex].x - viaBumps[overlapBumpIndex].x ; 
                ddy = offsetBumps[overlapBumpIndex].y - viaBumps[overlapBumpIndex].y ; 
                dx = tileNode2.crossedViaNode1->x - tileNode2.crossedViaNode2->x ; 
                dy = tileNode2.crossedViaNode1->y - tileNode2.crossedViaNode2->y ; 
                r1 = (tileNode2.crossedViaNode1->type==DUMMY) ? 0.0 : designRule.viaPadDiameter/2 ;
                r2 = (tileNode2.crossedViaNode1->type==DUMMY) ? 0.0 : designRule.viaPadDiameter/2 ;
                distance = sqrt(dx*dx + dy*dy) - sqrt(ddx*ddx + ddy*ddy) - r1 - r2 ; 
            }else{
                dx =  tileNode2.crossedViaNode1->x - tileNode2.crossedViaNode2->x ; 
                dy = tileNode2.crossedViaNode1->y - tileNode2.crossedViaNode2->y ; 
                r1 = (tileNode2.crossedViaNode1->type==DUMMY) ? 0.0 : designRule.viaPadDiameter/2 ;
                r2 = (tileNode2.crossedViaNode1->type==DUMMY) ? 0.0 : designRule.viaPadDiameter/2 ;
                distance = sqrt(dx*dx + dy*dy) - r1 - r2 ; 
            }

            
            *tileNode2.capacity = max(floor(distance/(designRule.minimumLineWidth + designRule.minimumLineSpacing)), 0.0) ; 
           
            // string name = "ground" ;
            // cout << tileNode2.crossedViaNode1->x << " " << tileNode2.crossedViaNode1->y << " " << tileNode2.crossedViaNode2->x << " " << tileNode2.crossedViaNode2->y << "\n" ; 
            // debugNets.push_back(Net(name, {{tileNode2.crossedViaNode1->x, tileNode2.crossedViaNode1->y, tileNode2.crossedViaNode2->x, tileNode2.crossedViaNode2->y}})) ;
            // dx = (tileNode2.crossedViaNode1->x + tileNode2.crossedViaNode2->x)/2 ;
            // dy = (tileNode2.crossedViaNode1->y + tileNode2.crossedViaNode2->y)/2 ;
            // debugLabels.push_back({to_string(*tileNode2.capacity), dx, dy}) ;
        }
    }
}
void Router::ConstructRoutingGraph(const vector<Bump>& routingBumps, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps, RoutingGraph2& graph, double minimumHorizontalSpace) {
    RoutingGraph2 RDL1, RDL2, RDL3 ; //RDL1 for left part, RDL2 for right part
    
    vector<Bump> die1Bumps, die2Bumps, die1PaddingBumps, die2PaddingBumps ; 
    vector<double> die1Coordinate, die2Coordinate ;  

    ofstream outputViaFile, outputTriangulationFile, outputNetlistFile ; 
    die1Coordinate = coordinate ; die1Coordinate[2] = die1Coordinate[2]/2 - minimumHorizontalSpace ; 
    die2Coordinate = coordinate ; die2Coordinate[0] = die2Coordinate[2]/2 + minimumHorizontalSpace ; 

    copy_if(routingBumps.begin(), routingBumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(offsetBumps.begin(), offsetBumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;

    copy_if(routingBumps.begin(), routingBumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;
    copy_if(offsetBumps.begin(), offsetBumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;

    PaddingBumps(die1Bumps, die1PaddingBumps, die1Coordinate, minimumHorizontalSpace) ;
    PaddingBumps(die2Bumps, die2PaddingBumps, die2Coordinate, minimumHorizontalSpace) ;

    // 建立Via nodes
    CreateViaNodes(RDL1, die1Bumps) ; CreateViaNodes(RDL1, die1PaddingBumps) ;
    CreateViaNodes(RDL2, die2Bumps) ; CreateViaNodes(RDL2, die2PaddingBumps) ;

    // 建立Via-via edges & Tile nodes & via-tile edges 使用CDT三角化法
    Triangulation(RDL1);
    Triangulation(RDL2);

    // 建立via-tile edges & tile-tile edges
    ConnectTileTileEdges(RDL1) ;     
    ConnectTileTileEdges(RDL2) ;     

    CombineRDLs(RDL1, RDL2, RDL3, die1Coordinate, die2Coordinate) ; 
    SetCapacity(RDL3, offsetBumps, viaBumps) ; 

    graph = RDL3 ;
}

void Router::GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph){

}


//------------------------------------------ Router Method End ------------------------------------------ 
