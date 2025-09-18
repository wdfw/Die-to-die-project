#include "GlobalRouting2.hpp"

//------------------------------------------ Router Method Begin ------------------------------------------ 

Router::Router(const DesignRule& designRule, const vector<Bump>& bumps, const vector<double>& coordinate){
    this->designRule = designRule ; 
    this->bumps = bumps ; 
    this->coordinate = coordinate ; 
}

void Router::Initial(){
    this->routingTimes.clear() ; 
    this->routingGraphs.clear() ; 
}

void Router::GlobalRoute(const string& outputDirectories){
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
            if(fabs(bump.y - minY) < epsilonY){
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
        ConstructRoutingGraph(directoryPath, layer, horizontalSpace) ;
        routingTimes.push_back(timer.GetDurationMilliseconds()); 

        leftBumpCount = 0 ; // for(const auto& bump : die1Bumps) if(bump.type==SIGNAL) ++leftBumpCount ;
        
    }
}


void Router::FindLeftmostInEachRow(const vector<Bump>& bumps, vector<Bump>& leftMostBumps) {
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
void Router::PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, const vector<double>& coordinate, double minimumHorizontalSpace){
    double min_x = coordinate[0], max_x = coordinate[2], current_x ;
    Bump dummyBump ;

    vector<Bump> leftmostBumps ; FindLeftmostInEachRow(bumps, leftmostBumps) ;
    for(auto bump:leftmostBumps) debugBumps.push_back(bump) ;
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
            set<shared_ptr<ViaNode2>> uniqueViaNodes ; 
            for(auto& viaNode : graph.tileNodes[i]->viaNodes) uniqueViaNodes.insert(viaNode) ; 
            for(auto& viaNode : graph.tileNodes[j]->viaNodes) uniqueViaNodes.insert(viaNode) ; 

            if(uniqueViaNodes.size()==4){

                graph.tileNodes[i]->tileNodes.push_back(graph.tileNodes[j]) ;
                graph.tileNodes[j]->tileNodes.push_back(graph.tileNodes[i]) ;

                graph.tileNodes[j]->tileNodes.back().capacity = graph.tileNodes[i]->tileNodes.back().capacity ; 

                // cout << graph.tileNodes[i]->tileNodes.back().capacity << " " << graph.tileNodes[j]->tileNodes.back().capacity << "\n"; 
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
void Router::CombineDie1Die2(RoutingGraph2& graph1, RoutingGraph2& graph2, const vector<double>& coordinate1, const vector<double>& coordinate2){
    double cx = (coordinate1[2] + coordinate2[0])/2 ;
    double cy = (coordinate1[1] + coordinate1[3] + coordinate2[1] + coordinate2[3])/4 ;
    shared_ptr<TileNode2> centralTileNode = make_shared<TileNode2>("Dummy", DUMMY, 0, cx, cy) ; 

    for(auto& graph : {graph1, graph2}){
        int n = graph.tileNodes.size(), m = graph.viaNodes.size() ; 
        for(int i=0; i<n; i++){
            int crossedCount = 0 ; 
            for(int j=0; j<m; j++){
                for(int k=0; k<graph.viaNodes[j]->viaNodes.size(); k++){
                    pair<double, double> pt1 = {centralTileNode->x, graph.tileNodes[i]->y} ;
                    pair<double, double> pt2 = {graph.tileNodes[i]->x, graph.tileNodes[i]->y} ;
                    pair<double, double> pt3 = {graph.viaNodes[j]->x, graph.viaNodes[j]->y} ;
                    pair<double, double> pt4 = {graph.viaNodes[j]->viaNodes[k]->x, graph.viaNodes[j]->viaNodes[k]->y} ;
                    double d1 = cross(pt1, pt2, pt3) ;
                    double d2 = cross(pt1, pt2, pt4) ; 
                    double d3 = cross(pt3, pt4, pt1) ;
                    double d4 = cross(pt3, pt4, pt2) ;
                    if( ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) or (d3 < 0 && d4 > 0)) ){
                        ++crossedCount ;
                    } 
                }
            }

            if(crossedCount<=2){
                graph.tileNodes[i]->tileNodes.push_back(centralTileNode) ;
                centralTileNode->tileNodes.push_back(graph.tileNodes[i]) ;
                centralTileNode->tileNodes.back().capacity = graph.tileNodes[i]->tileNodes.back().capacity ;
            }
        }
    }
}

void Router::ConstructRoutingGraph(const string& outputPath, int layer, double minimumHorizontalSpace){
    RoutingGraph2 RDL1, RDL2 ; //RDL1 for left part, RDL2 for right part
    
    vector<Bump> die1Bumps, die2Bumps, die1PaddingBumps, die2PaddingBumps ; 
    vector<double> die1Coordinate, die2Coordinate ;  

    ofstream outputViaFile, outputTriangulationFile, outputNetlistFile ; 
    die1Coordinate = coordinate ; die1Coordinate[2] = die1Coordinate[2]/2 - minimumHorizontalSpace ; 
    die2Coordinate = coordinate ; die2Coordinate[0] = die2Coordinate[2]/2 + minimumHorizontalSpace ; 

    copy_if(bumps.begin(), bumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(bumps.begin(), bumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;

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

    CombineDie1Die2(RDL1, RDL2, die1Coordinate, die2Coordinate) ; 
    // AddAccessViaEdges(RDL1) ;
    // AddAccessViaEdges(RDL2) ;

    // AddCrossTileEdges(RDL1, designRule) ;
    // AddCrossTileEdges2(RDL1, designRule) ;
    cout <<  outputPath + "via_layer_" + to_string(layer) << "\n" ;
    outputViaFile.open(outputPath + "via_layer_" + to_string(layer)) ;
    for(auto& bump : die1Bumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die2Bumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die1PaddingBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : die2PaddingBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    for(auto& bump : debugBumps) outputViaFile << Bump2Str(bump) << "\n" ;
    outputViaFile.close() ;

    
    outputTriangulationFile.open(outputPath + "triangulation_edge");

    for(auto& RDL : {RDL1, RDL2}){
        for(auto& start : RDL.viaNodes){
            for(auto& traget : start->viaNodes){
                outputTriangulationFile << start->x << " " << start->y << " " << traget->x << " " << traget->y << "\n" ;
            }
        }
        for(auto& start : RDL.tileNodes){
            for(auto& traget : start->tileNodes){
                outputTriangulationFile << start->x << " " << start->y << " " << traget->x << " " << traget->y << "\n" ;
            }
        }
        // for(auto& start : RDL.tileNodes){
        //     for(auto& traget : start->viaNodes){
        //         outputTriangulationFile << start->x << " " << start->y << " " << traget->x << " " << traget->y << "\n" ;
        //     }
        // }
    }
    // for(auto& viaNode : debugEdges) outputTriangulationFile << edgeNode.start.first << " " << edgeNode.start.second << " " << edgeNode.end.first << " " << edgeNode.end.second << "\n" ;
    outputTriangulationFile.close() ;

    // outputNetlistFile.open(outputPath + "netlist_" + to_string(layer));
    // for(auto& net : debugNets){
    //     for(auto& seg : net){
    //         outputNetlistFile << net.name << " " << get<0>(seg) << " " << get<1>(seg) << " " << get<2>(seg) << " " << get<3>(seg) << "\n" ;
    //     }
    // }
    // outputNetlistFile.close() ;
}

//------------------------------------------ Router Method End ------------------------------------------ 
