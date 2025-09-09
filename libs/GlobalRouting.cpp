#include "GlobalRouting.hpp"

const double EPSILON_Y = 1e-6 ;
const double EPSILON_X = 1e-3;
vector<Bump> debugBumps ; 
vector<EdgeNode> debugEdges ; 

int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routingCoordinate, 
                const DesignRule& designRule, const string& outputDirectories, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes){

    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> marginalBumps ;
    
    double horizontalSpace, minY = numeric_limits<double>::max() ; 
    
    string directoryPath ; 
    Timer timer;
    
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

    return globalRouteTimes.size() ; // return #layers
}

void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routingCoordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL){
    RoutingGraph RDL1, RDL2 ; 
    
    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> die1LeftmostBumps, die2LeftmostBumps ; // = FindLeftmostInEachRow(bumps);
    vector<Bump> die1DummyBumps, die2DummyBumps ; // = ProcessLeftAndRight(bumps, leftmostBumps, routingCoordinate, Hor_SPACING_X);
    vector<double> die1Coordinate, die2Coordinate ; 

    ofstream outputViaFile, outputTriangulationFile ; 
    //未來需要調整的部分!!!
    die1Coordinate = routingCoordinate ; die1Coordinate[2] = routingCoordinate[2]/2 - Hor_SPACING_X ; 
    die2Coordinate = routingCoordinate ; die2Coordinate[0] = routingCoordinate[2]/2 + Hor_SPACING_X ; 

    copy_if(bumps.begin(), bumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(bumps.begin(), bumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;

    die1LeftmostBumps = FindLeftmostInEachRow(die1Bumps) ;
    die2LeftmostBumps = FindLeftmostInEachRow(die2Bumps) ;

    die1DummyBumps = ProcessLeftAndRight(die1Bumps, die1LeftmostBumps, die1Coordinate, Hor_SPACING_X) ;
    die2DummyBumps = ProcessLeftAndRight(die2Bumps, die2LeftmostBumps, die2Coordinate, Hor_SPACING_X) ;

    for(auto bump : die1Bumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    for(auto bump : die1DummyBumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    for(auto bump : die2Bumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    for(auto bump : die2DummyBumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    // 使用 CDT 進行三角化
    RDL1.edge_nodes = Triangulation(RDL1.via_nodes, designRule);
    RDL2.edge_nodes = Triangulation(RDL2.via_nodes, designRule);
    
    
    FindAdjacentViaNodes(RDL1) ; //!!!!!
    FindAdjacentViaNodes(RDL2) ; //!!!!!

    AddAccessViaEdges(RDL1);
    AddAccessViaEdges(RDL2);

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
}

vector<EdgeNode*> FindAdjacentEdgeNodes(ViaNode &via, vector<EdgeNode> &edge_nodes) {
    vector<EdgeNode*> adjacent_edges;
    
    for (auto &edge : edge_nodes) {
        // 取出 edge 的兩端 ViaNode
        ViaNode &via1 = edge.vias[0];
        ViaNode &via2 = edge.vias[1];

        // **檢查 EdgeNode 是否與 ViaNode 相鄰**
        bool connects_to_via = ((via1.dieName == via.dieName && via1.type == via.type && via1.id == via.id) ||
                                (via2.dieName == via.dieName && via2.type == via.type && via2.id == via.id));
        if (connects_to_via) 
            adjacent_edges.push_back(&edge);
        
    }

    return adjacent_edges;
}

vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& coordinate, double Hor_SPACING_X){
    double min_x = coordinate[0], max_x = coordinate[2] ;

    Bump newPoint ;
    vector<Bump> dummies;

    double current_x ; 
    
    for (const auto& point : leftmostBumps) {
        double current_x = point.x;

        // 向左移動
        while (current_x - Hor_SPACING_X >= min_x + Hor_SPACING_X) {
            current_x -= Hor_SPACING_X ;
            newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
                dummies.push_back(newPoint);
            }
        }

        current_x = min_x ; 
        newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
            dummies.push_back(newPoint);
        }

        // 向右移動
        current_x = point.x;
        while (current_x + Hor_SPACING_X <= max_x - Hor_SPACING_X) {
            current_x += Hor_SPACING_X;
            newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
                dummies.push_back(newPoint);
            }
        }

        current_x = max_x ; 
        newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);
        if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
            dummies.push_back(newPoint);
        }
    }



    return dummies;
}


vector<Bump> FindLeftmostInEachRow(const vector<Bump>& bumps) {
    vector<vector<Bump>> y_groups;

    // 先依 y 排序，確保同 row 的點靠在一起
    vector<Bump> sorted_vias = bumps;
    sort(sorted_vias.begin(), sorted_vias.end(), [](const Bump& a, const Bump& b) {return a.y < b.y;});

    // EPSILON_Y 是 允許的 y 誤差範圍
    for (const auto& via : sorted_vias) {
        bool added = false;
        for (auto& group : y_groups) {
            if (!group.empty() && fabs(group[0].y - via.y) < EPSILON_Y) {
                group.push_back(via);
                added = true;
                break;
            }
        }
        if (!added) {
            y_groups.push_back({via});
        }
    }

    // 每組選出 x 最小的 bump（即最左點）
    vector<Bump> result;
    for (const auto& group : y_groups) {
        Bump leftmost = *min_element(group.begin(), group.end(), [](const Bump& a, const Bump& b) {return a.x < b.x;});
        result.push_back(leftmost);
    }

    return result;
}
bool f = false ; 
void FindAdjacentViaNodes(RoutingGraph &RDL){
    int n = RDL.edge_nodes.size() ; 

    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            for(int k=j+1; k<n; k++){
                
                map<pair<DieType, int>, int> countID, countY ; 
                map<pair<DieType, int>, ViaNode> uniqueVias ; 
                double aveX = 0.0, aveY = 0.0 ;
                bool isClosed = true ;
                
                ++ countID[{RDL.edge_nodes[i].vias[0].type, RDL.edge_nodes[i].vias[0].id}] ;
                ++ countID[{RDL.edge_nodes[j].vias[0].type, RDL.edge_nodes[j].vias[0].id}] ;
                ++ countID[{RDL.edge_nodes[k].vias[0].type, RDL.edge_nodes[k].vias[0].id}] ;

                ++ countID[{RDL.edge_nodes[i].vias[1].type, RDL.edge_nodes[i].vias[1].id}] ;
                ++ countID[{RDL.edge_nodes[j].vias[1].type, RDL.edge_nodes[j].vias[1].id}] ;
                ++ countID[{RDL.edge_nodes[k].vias[1].type, RDL.edge_nodes[k].vias[1].id}] ;

                uniqueVias[{RDL.edge_nodes[i].vias[0].type, RDL.edge_nodes[i].vias[0].id}] = RDL.edge_nodes[i].vias[0] ;
                uniqueVias[{RDL.edge_nodes[j].vias[0].type, RDL.edge_nodes[j].vias[0].id}] = RDL.edge_nodes[j].vias[0] ;
                uniqueVias[{RDL.edge_nodes[k].vias[0].type, RDL.edge_nodes[k].vias[0].id}] = RDL.edge_nodes[k].vias[0] ;

                uniqueVias[{RDL.edge_nodes[i].vias[1].type, RDL.edge_nodes[i].vias[1].id}] = RDL.edge_nodes[i].vias[1] ;
                uniqueVias[{RDL.edge_nodes[j].vias[1].type, RDL.edge_nodes[j].vias[1].id}] = RDL.edge_nodes[j].vias[1] ;
                uniqueVias[{RDL.edge_nodes[k].vias[1].type, RDL.edge_nodes[k].vias[1].id}] = RDL.edge_nodes[k].vias[1] ;

                aveX += RDL.edge_nodes[i].vias[0].x ; 
                aveX += RDL.edge_nodes[j].vias[0].x ; 
                aveX += RDL.edge_nodes[k].vias[0].x ; 

                aveX += RDL.edge_nodes[i].vias[1].x ; 
                aveX += RDL.edge_nodes[j].vias[1].x ; 
                aveX += RDL.edge_nodes[k].vias[1].x ; 

                aveY += RDL.edge_nodes[i].vias[0].y ; 
                aveY += RDL.edge_nodes[j].vias[0].y ; 
                aveY += RDL.edge_nodes[k].vias[0].y ;

                aveY += RDL.edge_nodes[i].vias[1].y ; 
                aveY += RDL.edge_nodes[j].vias[1].y ; 
                aveY += RDL.edge_nodes[k].vias[1].y ; 

                for(auto& [v, c] : countID){
                    if(c!=2) isClosed = false ; 
                }

                if(isClosed){   
                    // debugBumps.emplace_back("DUMMY", DUMMY, 0, aveX/6, aveY/6) ;
                    vector<ViaNode> vias ; 
                    for(auto& [c, v] : uniqueVias) vias.push_back(v) ;

                    vector<Edge> edges = {
                        {RDL.edge_nodes[i].id, RDL.edge_nodes[i].vias, RDL.edge_nodes[i].capacity},
                        {RDL.edge_nodes[j].id, RDL.edge_nodes[j].vias, RDL.edge_nodes[j].capacity},
                        {RDL.edge_nodes[k].id, RDL.edge_nodes[k].vias, RDL.edge_nodes[k].capacity},
                    } ;

                    RDL._edge_nodes.emplace_back(RDL._edge_nodes.size(), aveX/6, aveY/6, vias, edges) ;
                }
            }
        }
    }
}
void AddAccessViaEdges(RoutingGraph &RDL) {
    int current_access_via_id = static_cast<int>(RDL.access_via_edges.size());  // 記錄初始 ID
    set<tuple<string, DieType, int, int>> unique_access_via_edges; // 確保 (dieName, type, id, edge_id) 唯一


    for (auto &via : RDL.via_nodes) {
        vector<EdgeNode*> adjacent_edges = FindAdjacentEdgeNodes(via, RDL.edge_nodes);

        for (EdgeNode* edge : adjacent_edges) {
            tuple<string, DieType, int, int> access_via_tuple = {via.dieName, via.type, via.id, edge->id};

            // **如果已經存在，則跳過**
            if (unique_access_via_edges.find(access_via_tuple) != unique_access_via_edges.end()) {
                continue;
            }

            // **新增 AccessViaEdge**
            RDL.access_via_edges.emplace_back(current_access_via_id++, via, *edge);
            unique_access_via_edges.insert(access_via_tuple);
        }
    }

}

vector<EdgeNode*> findAdjacentEdgeNodes(ViaNode &via, vector<EdgeNode> &edge_nodes) {
    vector<EdgeNode*> adjacent_edges;
    
    for (auto &edge : edge_nodes) {
        // 取出 edge 的兩端 ViaNode
        ViaNode &via1 = edge.vias[0];
        ViaNode &via2 = edge.vias[1];

        // **檢查 EdgeNode 是否與 ViaNode 相鄰**
        bool connects_to_via = ((via1.dieName == via.dieName && via1.type == via.type && via1.id == via.id) ||
                                (via2.dieName == via.dieName && via2.type == via.type && via2.id == via.id));
        if (connects_to_via) 
            adjacent_edges.push_back(&edge);
        
    }

    return adjacent_edges;
}

vector<EdgeNode> Triangulation(const vector<ViaNode>& via_nodes, const DesignRule& designRule) {
    map<Point, ViaNode> point_to_via;
    vector<Point> allPoints;

    for (const auto &node : via_nodes) {
        Point p(node.x, node.y);
        point_to_via[p] = node;
        allPoints.push_back(p);
    }

    CDT cdt;

    // === Y 分組（允許誤差）===
    vector<vector<Point>> y_groups;
    sort(allPoints.begin(), allPoints.end(), [](const Point& a, const Point& b) {
        return a.y() < b.y();
    });

    for (const auto& pt : allPoints) {
        bool added = false;
        for (auto& group : y_groups) {
            if (!group.empty() && fabs(group[0].y() - pt.y()) < EPSILON_Y) {
                group.push_back(pt);
                added = true;
                break;
            }
        }
        if (!added) {
            y_groups.push_back({pt});
        }
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
    vector<EdgeNode> edges;
    int edgeId = 0;
    
    // ofstream outFile(output_path + "triangulation_edge", ofstream::app);

    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        auto segment = cdt.segment(*edge);
        double dx = segment.source().x() - segment.target().x();
        double dy = segment.source().y() - segment.target().y();
        double L = sqrt(dx * dx + dy * dy);
        double ux = dx / L, uy = dy / L;

        double newSourceX, newSourceY, newTargetX, newTargetY;
        if (point_to_via[segment.source()].type == DUMMY) {
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

        if (point_to_via[segment.target()].type == DUMMY) {
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
        int capacity = floor(distance / (2*(designRule.minimumLineWidth + designRule.minimumLineSpacing))); 
        if(capacity > 5 || capacity < 5) capacity = 5; // 限制最大容量為 5

        EdgeNode newEdgeNode(edgeId++, x, y,
            vector<ViaNode>{point_to_via[segment.source()], point_to_via[segment.target()]}, // Edge Node 紀錄連接的兩個 ViaNode
            capacity);

        string src_type = point_to_via[segment.source()].type + "_" + to_string(point_to_via[segment.source()].id);
        string tgt_type = point_to_via[segment.target()].type + "_" + to_string(point_to_via[segment.target()].id);

        // 判斷 net_sequence 的方向
        if (fabs(segment.source().y() - segment.target().y()) > EPSILON_Y) { // Edge Node 連接的兩個 ViaNode 的 y 值不同
            if (segment.source().y() < segment.target().y()) {
                newEdgeNode.net_sequence = {src_type, tgt_type}; // 紀錄 net sequence list 的 head node 和 tail node, 並同時紀錄兩端的座標, detailed route 較方便
                newEdgeNode.start = {newSourceX, newSourceY};
                newEdgeNode.end   = {newTargetX, newTargetY};
            } else {  // source.y > target.y
                newEdgeNode.net_sequence = {tgt_type, src_type};
                newEdgeNode.start = {newTargetX, newTargetY};
                newEdgeNode.end   = {newSourceX, newSourceY};
            }
        } else {
            // Edge Node 連接的兩個 ViaNode 的 y 值相同, 比較 x 值
            if (segment.source().x() > segment.target().x()) {
                newEdgeNode.net_sequence = {src_type, tgt_type};
                newEdgeNode.start = {newSourceX, newSourceY};
                newEdgeNode.end   = {newTargetX, newTargetY};
            } else {
                newEdgeNode.net_sequence = {tgt_type, src_type};
                newEdgeNode.start = {newTargetX, newTargetY};
                newEdgeNode.end   = {newSourceX, newSourceY};
            }
        }

        edges.push_back(newEdgeNode);
    }

    // outFile.close();
    return edges;
}
