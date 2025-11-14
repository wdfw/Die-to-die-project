#include "Router.hpp"

vector<Bump> debugBumps ; 
// vector<EdgeNode> debugEdges ; 
vector<Net> debugNets ; 
vector<Teardrop> debugTeardrops ; 
vector<tuple<string, double, double>> debugLabels ; 
vector<Timer> globalTimers(100) ; 
vector<double> globalExecTimes(100, 0.0) ; 
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
void Router::SelectRoutingBumps(const vector<Bump>& bumps, vector<Bump>& routingBumps, vector<Bump>& offsetBumps, int selectNum){
    //Die1 越靠近中線的越先繞, feedback = 預估可繞線的數量

    int leftSignalCount = 0 ;
    set<int> selectedID ; 
    vector<int> allSignalBumpIndexs ;

    // else if(feedback==2.0) selectedID = set<int>{27, 28, 29, 30, 31, 32, 33, 34} ;
    // else if(feedback==3.0) selectedID = set<int>{19, 20, 21, 22, 23, 24, 25, 26} ; 
    // selectedID = set<int>{40, 36, 41, 37,42,34,35} ;

    for(int i=0; i<bumps.size(); ++i){
        if(bumps[i].type==SIGNAL && bumps[i].name=="DIE1") allSignalBumpIndexs.push_back(i) ; 
    }
    
    sort(allSignalBumpIndexs.begin(), allSignalBumpIndexs.end(), [&bumps](int i, int j){return bumps[i].x>bumps[j].x;} ) ;
    for(int i=0; i<min(selectNum, int(allSignalBumpIndexs.size())); ++i) selectedID.insert(bumps[allSignalBumpIndexs[i]].id) ;
    routingBumps.clear() ; 
    offsetBumps.clear() ; 
    for(int i=0; i<bumps.size(); i++){
        if(bumps[i].type==VSS){
            routingBumps.push_back(bumps[i]) ;
            offsetBumps.push_back(bumps[i]) ;
        }else if(bumps[i].type==VDD){
            routingBumps.push_back(bumps[i]) ;
            offsetBumps.push_back(bumps[i]) ; //!!!!!!!!
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
    vector<double> currentCorrdinate = coordinate ; 
    vector<GlobalNet> globalNets ;
    vector<DetailedNet> detailedNets ;

    for(int layer = 1, leftBumpCount = bumps.size(); unroutedBumps.size() ; ++layer){
        RoutingGraph2 routingGraph ; 
        directoryPath = outputDirectories + "RDL" + to_string(layer) + "/";
        int selectedNum = 0 ; 
        double feedback = 0.0 ;

         for(int i=0; i<unroutedBumps.size(); ++i){
            if(unroutedBumps[i].type==SIGNAL && unroutedBumps[i].name=="DIE1"){
                ++selectedNum ; 
            }
        }
        selectedNum = min(selectedNum, 10) ; //!!!!!!
       
        if (!filesystem::exists(directoryPath)) filesystem::create_directories(directoryPath);
        
        timer.SetClock() ;
        FindHorizontalSpace(unroutedBumps, horizontalSpace) ;
        
        cout << "----------------------------------Layer: " << layer << "-----------------------------------\n" ; 

        do{
            debugNets.clear() ; debugLabels.clear() ; debugTeardrops.clear() ; 
            SelectRoutingBumps(unroutedBumps, routingBumps, offsetBumps, selectedNum) ;
            CreateViaBumps(offsetBumps, viaBumps) ;
            ConstructRoutingGraph(routingBumps, offsetBumps, viaBumps, routingGraph, horizontalSpace, currentCorrdinate) ;
            
             
            // for(auto& tileNode : routingGraph.tileNodes) debugLabels.push_back(tuple<string, double, double>{to_string(tileNode->id), tileNode->x, tileNode->y}) ; 
            
             
            feedback = GlobalRoute(routingBumps, routingGraph, globalNets) ; 
            selectedNum += feedback ;
        }while(feedback!=0) ; 

        DetailRoute(routingBumps, routingGraph, globalNets, detailedNets) ; 

        GenerateGraphFile(routingBumps, offsetBumps, viaBumps, routingGraph, layer, directoryPath) ;
        unroutedBumps = viaBumps ; 

        currentCorrdinate[0] -=  designRule.minimumViaSpacing + designRule.viaOpeningDiameter ; 
        currentCorrdinate[2] -=  designRule.minimumViaSpacing + designRule.viaOpeningDiameter ; 
        routingTimes.push_back(timer.GetDurationMilliseconds()); 
    }
}

void Router::CreateViaBumps(const vector<Bump>& offsetBumps, vector<Bump>& viaBumps){
    viaBumps.clear() ; 
    for(auto& bump : offsetBumps){
        viaBumps.push_back(bump) ; 
        viaBumps.back().x = bump.x - ( designRule.minimumViaSpacing + designRule.viaOpeningDiameter) ;
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

    outputNetlistFile.open(directoryPath + "debug_net");
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


    
    outputDebugLabelFile.open(directoryPath + "teardrop_" + to_string(layer)) ;
    for(auto& debugTeardrop : debugTeardrops){
        outputDebugLabelFile << debugTeardrop.name << " " << DieType2Str(debugTeardrop.type) << " " << debugTeardrop.id << " " << debugTeardrop.sx << " " << debugTeardrop.sy << " " << debugTeardrop.ex << " " << debugTeardrop.ey << "\n" ;
    }
    outputDebugLabelFile.close() ;
}

void Router::PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, const vector<double>& coordinate, double minimumHorizontalSpace){
    vector<Bump> bumpStacks = bumps ;
    vector<Bump> edgeBumps ;
    bool isOnEdge, isSearched ;
    double r = minimumHorizontalSpace ; 

    double cos60 = 1.0/2 ; 
    double sin60 = sqrt(3)/2 ; 
    for(int i=0; i<bumpStacks.size(); ++i){
        const Bump& bump = bumpStacks[i] ;
        double cx = bump.x, cy = bump.y ; 
        for(auto [rx, ry] : vector<pair<double,double>>{{1,0},{-1,0},{cos60,sin60},{-cos60,sin60},{cos60,-sin60},{-cos60,-sin60}}){
            double nx = cx+rx*r, ny = cy+ry*r ; 

            if(ny<coordinate[1] || ny>coordinate[3]) break ;

            if(nx<coordinate[0]){
                nx = coordinate[0] ; isOnEdge = true ; 
            }else if(nx>coordinate[2]){
                nx = coordinate[2] ; isOnEdge = true ; 
            }else{
                isOnEdge = false ; 
            }

            isSearched = false ; 

            
            for(auto& serchedBump : bumps){
                if(serchedBump==bump) continue ;

                if(abs(serchedBump.x-nx)<0.5 && abs(serchedBump.y-ny)<0.5){
                    isSearched = true ; break;
                }
            }

            if(!isSearched){
                for(auto& serchedBump : paddingBumps){
                    if(serchedBump==bump) continue ;

                    if(abs(serchedBump.x-nx)<0.5 && abs(serchedBump.y-ny)<0.5){
                        isSearched = true ; break;
                    }
                }
            }

            if(!isSearched){
                paddingBumps.push_back(Bump("Dummy", DUMMY, paddingBumps.size(), nx, ny)) ;
                if(isOnEdge) edgeBumps.push_back(paddingBumps.back()) ; 
                else bumpStacks.push_back(paddingBumps.back()) ; 
            }            
        }
    }

    for(int i=0; i<edgeBumps.size(); ++i){
        for(int j=0; j<paddingBumps.size(); ++j){
            if(edgeBumps[i]==paddingBumps[j]) continue ;
            if(abs(paddingBumps[j].x-edgeBumps[i].x)<0.5*r && abs(paddingBumps[j].y-edgeBumps[i].y)<0.5){
                swap(paddingBumps[j], paddingBumps.back()) ; 
                paddingBumps.pop_back() ; 
                --j ; 
            }
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
                

                if( crossedViaNodes[0]->name == crossedViaNodes[1]->name &&
                    crossedViaNodes[0]->id == crossedViaNodes[1]->id &&
                    crossedViaNodes[0]->type == crossedViaNodes[1]->type){
                        continue; 
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

                dieEscapeTileNodes.back()->tileNodes.back().crossedViaNode1 = (crossedViaNodes[0]->y < crossedViaNodes[1]->y) ? crossedViaNodes[0] : crossedViaNodes[1] ;
                dieEscapeTileNodes.back()->tileNodes.back().crossedViaNode2 = (crossedViaNodes[0]->y < crossedViaNodes[1]->y) ? crossedViaNodes[1] : crossedViaNodes[0] ;
            }
        }
    }

    sort(die1EscapeTileNodes.begin(), die1EscapeTileNodes.end(), [](const shared_ptr<TileNode2>& a, const shared_ptr<TileNode2>& b){return a->y<b->y; }) ;
    sort(die2EscapeTileNodes.begin(), die2EscapeTileNodes.end(), [](const shared_ptr<TileNode2>& a, const shared_ptr<TileNode2>& b){return a->y<b->y; }) ;

    if(die1EscapeTileNodes.size()!=die2EscapeTileNodes.size()) throw runtime_error("Invaild size on dies escape tile nodes (" + to_string(die1EscapeTileNodes.size()) + "," + to_string(die2EscapeTileNodes.size()) + ")") ; 

    for(int i=0; i<die1EscapeTileNodes.size(); ++i){
        if(i==0) die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1->viaNodes.push_back(die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1) ;
        die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2->viaNodes.push_back(die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2) ;
    }

    shared_ptr<ViaNode2> crossedViaNode1, crossedViaNode2 ; 
    vector<shared_ptr<TileNode2>> intermediateTileNodes ;
    for(int i=0; i<die1EscapeTileNodes.size(); i++){
        shared_ptr<TileNode2> intermediateTileNode = make_shared<TileNode2>("DUMMY", DUMMY, graph1.tileNodes.size()+graph2.tileNodes.size()+intermediateTileNodes.size(), 
                                                    (die1EscapeTileNodes[i]->x + die2EscapeTileNodes[i]->x)/2, (die1EscapeTileNodes[i]->y + die2EscapeTileNodes[i]->y)/2) ; 

        intermediateTileNode->tileNodes.push_back(die1EscapeTileNodes[i]) ;
        intermediateTileNode->tileNodes.back().crossedViaNode1 = die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 ; 
        intermediateTileNode->tileNodes.back().crossedViaNode2 = die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 ; 
        intermediateTileNode->tileNodes.back().capacity = die1EscapeTileNodes[i]->tileNodes.back().capacity ;

        die1EscapeTileNodes[i]->tileNodes.back() = intermediateTileNode ;
        die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 = intermediateTileNode->tileNodes.back().crossedViaNode1 ;
        die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 = intermediateTileNode->tileNodes.back().crossedViaNode2 ;
        die1EscapeTileNodes[i]->tileNodes.back().capacity = intermediateTileNode->tileNodes.back().capacity ;

        // cout << die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 << " " << die1EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 << "\n" ; 

        intermediateTileNode->tileNodes.push_back(die2EscapeTileNodes[i]) ;
        intermediateTileNode->tileNodes.back().crossedViaNode1 = die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 ; 
        intermediateTileNode->tileNodes.back().crossedViaNode2 = die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 ; 
        intermediateTileNode->tileNodes.back().capacity = die2EscapeTileNodes[i]->tileNodes.back().capacity ;

        die2EscapeTileNodes[i]->tileNodes.back() = intermediateTileNode ;
        die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode1 = intermediateTileNode->tileNodes.back().crossedViaNode1 ;
        die2EscapeTileNodes[i]->tileNodes.back().crossedViaNode2 = intermediateTileNode->tileNodes.back().crossedViaNode2 ;
        die2EscapeTileNodes[i]->tileNodes.back().capacity = intermediateTileNode->tileNodes.back().capacity ;

        intermediateTileNodes.push_back(intermediateTileNode) ; 
    }


    graph3 = RoutingGraph2(graph1) ; 
    for(auto& viaNode : graph2.viaNodes) graph3.viaNodes.push_back(viaNode) ;
    for(auto& tileNode : graph2.tileNodes) graph3.tileNodes.push_back(tileNode) ;
    for(auto& tileNode : intermediateTileNodes) graph3.tileNodes.push_back(tileNode) ;
}

void Router::SetCapacity(RoutingGraph2& graph){
    double dx, dy, ddx, ddy, distance, r1, r2 ;
    int overlapBumpIndex ; 

    for(auto& viaNode : graph.viaNodes){
        for(auto& tileNode : viaNode->tileNodes){
            (*tileNode.capacity) = numeric_limits<int>::max();
        }
    }

    map<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>, shared_ptr<EdgeNode>> edgeNodeMapping ; 
    for(auto& edgeNode : graph.edgeNodes) edgeNodeMapping[{edgeNode->viaNode1, edgeNode->viaNode2}] = edgeNode ; 

    for(auto& tileNode1 : graph.tileNodes){
        for(auto& tileNode2 : tileNode1->tileNodes){
            pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>> key1 = {tileNode2.crossedViaNode1, tileNode2.crossedViaNode2} ;
            pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>> key2 = {tileNode2.crossedViaNode2, tileNode2.crossedViaNode1} ;
            int capacity = (edgeNodeMapping.find(key1)!=edgeNodeMapping.end()) ? edgeNodeMapping[key1]->positions.size() : edgeNodeMapping[key2]->positions.size() ; 
            (*tileNode2.capacity) = max((capacity-1)/2, 0) ;
        }
    }
}
void Router::ConstructRoutingGraph(const vector<Bump>& routingBumps, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps, RoutingGraph2& graph, 
                                   double minimumHorizontalSpace, const vector<double>& coordinate) {
    RoutingGraph2 RDL1, RDL2, RDL3 ; //RDL1 for left part, RDL2 for right part
    
    vector<Bump> die1Bumps, die2Bumps, die1PaddingBumps, die2PaddingBumps ; 
    vector<double> die1Coordinate, die2Coordinate ;  

    ofstream outputViaFile, outputTriangulationFile, outputNetlistFile ; 
    die1Coordinate = coordinate ; die1Coordinate[2] = (die1Coordinate[2]+die1Coordinate[0])/2 - minimumHorizontalSpace ; 
    die2Coordinate = coordinate ; die2Coordinate[0] = (die2Coordinate[2]+die2Coordinate[0])/2 + minimumHorizontalSpace ; 

    copy_if(routingBumps.begin(), routingBumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(offsetBumps.begin(), offsetBumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;

    copy_if(routingBumps.begin(), routingBumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;
    copy_if(offsetBumps.begin(), offsetBumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;


    PaddingBumps(die1Bumps, die1PaddingBumps, die1Coordinate, minimumHorizontalSpace) ;
    PaddingBumps(die2Bumps, die2PaddingBumps, die2Coordinate, minimumHorizontalSpace) ;

    copy_if(viaBumps.begin(), viaBumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(viaBumps.begin(), viaBumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;


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
    CreateEdgeNodes(RDL3, offsetBumps, viaBumps) ; 
    SetCapacity(RDL3) ; 
    graph = RDL3 ;
}

void Router::CreateEdgeNodes(RoutingGraph2& graph, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps){

    shared_ptr<ViaNode2> viaNode1, viaNode2 ; 
    vector<pair<double, double>> positions ; 
    EdgeType type ; 
    double rx, ry ; 
    double bx, by ; 
    double dx, dy ; 
    double theta ; 
    double distance ;
    int capacity ;
    int offsetIndex ;
    set<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>> usedPairs ; 
    for(auto& tileNode : graph.tileNodes){
        for(auto& adjacentTileNode : tileNode->tileNodes){
            pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>> ptrPair = {adjacentTileNode.crossedViaNode1 , adjacentTileNode.crossedViaNode2} ;
            int capacity = *adjacentTileNode.capacity ; 
            if(usedPairs.find(ptrPair)==usedPairs.end()){
                usedPairs.insert(ptrPair) ;
                
                if(fabs(ptrPair.first->y - ptrPair.second->y) < 0.5) type = BaseEdge ; 
                else type = LegEdge ; 

                if(type==BaseEdge && ptrPair.first->x > ptrPair.second->x) swap(ptrPair.first, ptrPair.second) ; 
                else if(type==LegEdge && ptrPair.first->y > ptrPair.second->y) swap(ptrPair.first, ptrPair.second) ; 
                

                dx = ptrPair.second->x - ptrPair.first->x ; 
                dy = ptrPair.second->y - ptrPair.first->y ; 
                bx = ptrPair.first->x ;
                by = ptrPair.first->y ;

                // if(type==BaseEdge){ //處理offset via
                //     for(int i=0; i<offsetBumps.size(); ++i){
                //         auto& bump = offsetBumps[i] ;
                //         if(bump.name==ptrPair.second->name && bump.id==ptrPair.second->id && bump.type==ptrPair.second->type){
                //             dx = viaBumps[i].x - ptrPair.first->x ; dy = viaBumps[i].y - ptrPair.first->y ; 
                //             break;
                //         }
                //     }
                // }
                theta = atan2(dy, (dx)+1e-6) ; 
                

                bx += designRule.minimumLineWidth*cos(theta) ; dx -= 2*designRule.minimumLineWidth*cos(theta) ; 
                by += designRule.minimumLineWidth*sin(theta) ; dy -= 2*designRule.minimumLineWidth*sin(theta) ;

                if(ptrPair.first->type!=DUMMY){
                    bx += designRule.minimumViaPadSpacing/2*cos(theta) ; dx -= designRule.minimumViaPadSpacing/2*cos(theta) ; 
                    by += designRule.minimumViaPadSpacing/2*sin(theta) ; dy -= designRule.minimumViaPadSpacing/2*sin(theta) ;
                }
                
                if(ptrPair.second->type!=DUMMY){
                    dx -= designRule.minimumViaPadSpacing/2*cos(theta)  ; 
                    dy -= designRule.minimumViaPadSpacing/2*sin(theta) ;
                }

                if(type==BaseEdge && dx<=0) distance = 0 ;
                else if(type==LegEdge && dy<=0) distance = 0 ;
                else distance = sqrt(dx*dx+dy*dy) ; 
               
                capacity = max(floor(distance/(designRule.minimumLineWidth + designRule.minimumLineSpacing)), 0.0) ; 

                if(capacity!=0){
                    rx = dx/capacity ; ry = dy/capacity ; 
                    // cout << rx << " " << (designRule.minimumLineWidth + designRule.minimumLineSpacing)*cos(theta) << "\n" ;
                    bx += 0.5*rx ; by += 0.5*ry ;
                    // debugNets.push_back(Net(to_string(100000+debugNets.size()), {{bx, by, dx+bx, dy+by}})) ;
                }
                debugLabels.push_back({to_string(capacity), dx/2+bx, dy/2+by}) ;
                positions.clear() ; 
                
                
                for(int i=0; i<capacity; ++i){
                    positions.push_back({bx+rx*i, by+ry*i}) ; 
                }

                if(capacity!=0){
                    // debugNets.push_back(Net(to_string(100000+debugNets.size()), {{positions[0].first, positions[0].second, positions.back().first, positions.back().second}})) ;
                }

                graph.edgeNodes.push_back(make_shared<EdgeNode>(ptrPair.first, ptrPair.second, type, positions)) ;
            }
        }
    }
}
double Router::GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph, vector<GlobalNet>& nets){
    return 0 ; 
}

double Router::DetailRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph, vector<GlobalNet>& globalNets, vector<DetailedNet>& detailedNets){
    detailedNets.clear() ; 
    
    map<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>, shared_ptr<EdgeNode>> viaNodeToEdgeNode ;  

    map<TileToTileEdge, shared_ptr<EdgeNode>> edgeNodeMapping ;  
    map<shared_ptr<EdgeNode>, int> usedNum ;
    map<shared_ptr<EdgeNode>, DieType> lastUsedType ;

    for(auto& edgeNode : graph.edgeNodes){
        viaNodeToEdgeNode[{edgeNode->viaNode1, edgeNode->viaNode2}] = edgeNode ; 
        viaNodeToEdgeNode[{edgeNode->viaNode2, edgeNode->viaNode1}] = edgeNode ; 
    }

    for(auto& tileNode : graph.tileNodes){
        // cout << *tileNode->viaNodes[0] << " " << *tileNode->viaNodes[1] << " " << *tileNode->viaNodes[2] << "|\n" ; 
        for(auto& tileToTileEdge : tileNode->tileNodes){
            edgeNodeMapping[tileToTileEdge] = viaNodeToEdgeNode[{tileToTileEdge.crossedViaNode1,tileToTileEdge.crossedViaNode2}] ;
            // cout << *tileToTileEdge.crossedViaNode1 << " " << *tileToTileEdge.crossedViaNode2 << "\n" ; 
            // auto& vv = viaNodeToEdgeNode[{tileToTileEdge.crossedViaNode1,tileToTileEdge.crossedViaNode2}] ;
            // cout << *vv->viaNode1 << " " << *vv->viaNode2 << "|\n" ;
        }
        // cout << "\n\n" ;
    }

    for(int i=0; i<globalNets.size(); ++i){
        GlobalNet& gnet = globalNets[i] ;
        DetailedNet dnet ; 
        
        dnet.name = gnet.name ; 
        dnet.type = gnet.type ; 
        dnet.id = gnet.id ; 
        dnet.startViaNode = gnet.startViaNode ;
        dnet.endViaNode = gnet.endViaNode ;
        
        // cout << gnet.name << " " << DieType2Str(gnet.type) << " " << gnet.id << "\n" ;

        for(int j=1; j<gnet.size(); ++j){ //1開始是因為第0個TileToTile Node是空2邊的
            shared_ptr<EdgeNode> edgeNode = viaNodeToEdgeNode[{gnet[j].crossedViaNode1, gnet[j].crossedViaNode2}] ;

            if(lastUsedType[edgeNode]==VSS && !gnet.startViaNode) --usedNum[edgeNode] ; 

            if(usedNum.find(edgeNode)==usedNum.end()) usedNum[edgeNode] = 0 ;

            if(edgeNode->type==LegEdge) dnet.push_back({edgeNode, usedNum[edgeNode]}) ; 
            else{
                
                if(edgeNode->positions[0].second < dnet.back().first->positions[0].second) dnet.push_back({edgeNode, usedNum[edgeNode]}) ; //斜上
                else dnet.push_back({edgeNode, edgeNode->positions.size()-1-usedNum[edgeNode]}) ; // 斜下
            }
            
            ++usedNum[edgeNode] ; 

            if(gnet.startViaNode) lastUsedType[edgeNode] = SIGNAL ; 
            else lastUsedType[edgeNode] = VSS ; 
        }
        // cout << "\n" ;
        detailedNets.push_back(dnet) ; 
    }

    vector<string> colorTypes = {"blue", "red", "gray", "white"} ; 
        

    for(int i=0, v=0; i<detailedNets.size(); ++i, ++v){
        string ccolor = colorTypes[v%colorTypes.size()] ; 
        auto& detailedNet = detailedNets[i] ; 

        if(!detailedNet.startViaNode) ccolor = "green" ;

        if(detailedNet.startViaNode){
            debugNets.push_back(Net(ccolor, {{detailedNet.startViaNode->x , detailedNet.startViaNode->y,
                detailedNet[0].first->positions[detailedNet[0].second].first,  detailedNet[0].first->positions[detailedNet[0].second].second}})) ;
            
            double dx = detailedNet[0].first->positions[detailedNet[0].second].first-detailedNet.startViaNode->x ;
            double dy = detailedNet[0].first->positions[detailedNet[0].second].second-detailedNet.startViaNode->y ;
            double theta = atan2(dy, dx+1e-6) ;

            debugTeardrops.emplace_back( 
                detailedNet.startViaNode->name, 
                detailedNet.startViaNode->type,
                detailedNet.startViaNode->id,
                detailedNet.startViaNode->x,
                detailedNet.startViaNode->y,
                detailedNet.startViaNode->x + designRule.minimumTeardropDist*cos(theta),
                detailedNet.startViaNode->y + designRule.minimumTeardropDist*sin(theta)
            ) ;
        }

        // cout << detailedNets[i].size() << "\n" ;
        for(int j=0; j<detailedNets[i].size()-1; ++j){
            // cout << detailedNet[i].first->positions[detailedNet[i].second].first << " " <<  detailedNet[i].first->positions[detailedNet[i].second].second << "\n" ;
            debugNets.push_back(Net(ccolor, {{detailedNet[j].first->positions[detailedNet[j].second].first,  detailedNet[j].first->positions[detailedNet[j].second].second,
                    detailedNet[j+1].first->positions[detailedNet[j+1].second].first,  detailedNet[j+1].first->positions[detailedNet[j+1].second].second}})) ;
        }
        // cout << "\n" ;
        
        
        if(detailedNet.endViaNode){
            debugNets.push_back(Net(ccolor, {{detailedNet.back().first->positions[detailedNet.back().second].first,  detailedNet.back().first->positions[detailedNet.back().second].second,
                    detailedNet.endViaNode->x , detailedNet.endViaNode->y}})) ;

                      
            double dx =  detailedNet.endViaNode->x - detailedNet.back().first->positions[detailedNet.back().second].first ;
            double dy =  detailedNet.endViaNode->y - detailedNet.back().first->positions[detailedNet.back().second].second ; 
            double theta = atan2(dy, dx+1e-6) ;

            debugTeardrops.emplace_back( 
                detailedNet.endViaNode->name, 
                detailedNet.endViaNode->type,
                detailedNet.endViaNode->id,
                detailedNet.endViaNode->x,
                detailedNet.endViaNode->y,
                detailedNet.endViaNode->x - designRule.minimumTeardropDist*cos(theta),
                detailedNet.endViaNode->y - designRule.minimumTeardropDist*sin(theta)
            ) ;
        }
    }



    return 0 ; 
}

//------------------------------------------ Router Method End ------------------------------------------ 
