#include "GlobalRouting.hpp"

const double EPSILON_Y = 1e-6 ;
const double EPSILON_X = 1e-3;
vector<Bump> debugBumps ; 
vector<EdgeNode> debugEdges ; 
vector<Net> debugNets ; 
vector<tuple<string, double, double>> debugLabels ; 

int coutFlag = 0 ;
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


pair<vector<vector<ViaNode>>, vector<vector<ViaNode>>> findUpperLowerOtherPerRow(const vector<ViaNode>& vias, unordered_map<string, unordered_set<string>>& net_group_map) {
    map<double, vector<ViaNode>> row_nodes;  // 以 y 座標分類
    for (const auto& via : vias) {
        if (via.dieName == "DIE1" && via.type == SIGNAL)
            row_nodes[via.y].push_back(via);
    }

    // **結果：分成兩部分**
    vector<vector<ViaNode>> rightmost_list;  // 最右邊的兩個
    vector<vector<ViaNode>> remaining_list;  // 其餘的 ViaNode

    for (auto& [y, nodes] : row_nodes) {
        // **按 x 值排序**
        sort(nodes.begin(), nodes.end(), [](const ViaNode& a, const ViaNode& b) {
            return a.x > b.x;  // 由右向左排序
        });

        int num_nodes = nodes.size();
        unordered_set<string> row_net_set;  // **存這一列所有 `net_xx`，避免重複處理**
        
        for (const auto& node : nodes) {
            string net_name = "Net_" + to_string(node.id);
            row_net_set.insert(net_name);
        }

        // **將這列的 net_xx 加入同一組**
        for (const auto& net : row_net_set) {
            net_group_map[net] = row_net_set;  // **將整組直接存入**
        }

        if (num_nodes > 2) {
            // **將最右邊兩個加入 rightmost_list**
            rightmost_list.push_back({nodes[0], nodes[1]});

            // **其餘的加入 remaining_list**
            remaining_list.push_back(vector<ViaNode>(nodes.begin() + 2, nodes.end()));
        } else {
            // **如果只有兩個或更少，全部歸入 rightmost_list**
            rightmost_list.push_back(nodes);
            remaining_list.push_back({});
        }
    }

    return {rightmost_list, remaining_list};
}


// 根據 direction 找到離 via 最近的 EdgeNode
EdgeNode* findNearestEdgeNode(const ViaNode& via, vector<EdgeNode>& edge_nodes, const string& direction) { 
    double min_distance = numeric_limits<double>::max();
    EdgeNode* nearest_edge = nullptr;

    for (auto& edge : edge_nodes) {
        // **確保 EdgeNode 其中一端是該 ViaNode**
        bool has_via = (edge.vias[0].id == via.id && edge.vias[0].type == via.type) ||
                       (edge.vias[1].id == via.id && edge.vias[1].type == via.type);

        if (!has_via) continue;  // **如果 EdgeNode 兩端都不是該 ViaNode，則跳過**

        // **計算 EdgeNode 的中心點**
        double edge_x = (edge.vias[0].x + edge.vias[1].x) / 2.0;
        double edge_y = (edge.vias[0].y + edge.vias[1].y) / 2.0;

        // **計算與 ViaNode 的距離**
        double distance = sqrt(pow(edge_x - via.x, 2) + pow(edge_y - via.y, 2));

        // **篩選符合方向條件的 EdgeNode**
        bool valid_edge = false;
        if (direction == "top-right" && edge_x > via.x && edge_y < via.y) {
            valid_edge = true;
        } else if (direction == "bottom-right" && edge_x > via.x && edge_y > via.y) {
            valid_edge = true;
        } else if (direction == "bottom-left" && edge_x < via.x && edge_y > via.y) {
            valid_edge = true;
        }


        if (valid_edge && distance < min_distance) {
            min_distance = distance;
            nearest_edge = &edge;
        }
    }

    return nearest_edge;
}

CrossTileEdge* findLargestAngleCrossTileEdge(EdgeNode* current_edge, vector<CrossTileEdge>& cross_tile_edges, RoutingGraph& RDL) {
    CrossTileEdge* best_edge = nullptr;
    double max_angle = numeric_limits<double>::min();
    double max_x = numeric_limits<double>::min();
    const double EPSILON = 1e-6;

    for (auto& edge : cross_tile_edges) {
        if (edge.edges[0].id == current_edge->id || edge.edges[1].id == current_edge->id) {
            // **取得另一個 EdgeNode**
            EdgeNode* next_edge = (edge.edges[0].id == current_edge->id) ? &edge.edges[1] : &edge.edges[0];
            
            double next_x = next_edge->x;  // 找出 X 最大的值
            
            //  優先比較 X 最大
            if (next_x > max_x + EPSILON) {
                max_x = next_x;
                max_angle = edge.angle;
                best_edge = &edge;
            } else if (fabs(next_x - max_x) < EPSILON && edge.angle > max_angle) { //  若 X 相同（在誤差內），再選擇角度最大 (目前用不到, 因為所有角度都相同)
                max_angle = edge.angle;
                best_edge = &edge;
            }
        }
    }

    // cout<< "max_angle: "<<max_angle<<endl;

    return best_edge;
}

void routeUpperLower(vector<vector<ViaNode>>& rightmost_per_row, RoutingGraph& RDL) {
    for (size_t row = 0; row < 1; row++) {  //rightmost_per_row.size()
        int num_nodes = rightmost_per_row[row].size();
        if (num_nodes == 0) continue; 

        for (size_t i = 0; i < 1; i++) { //num_nodes
            ViaNode via = rightmost_per_row[row][i];

            cout << "Routing ViaNode: (Row: " << row << ", ID: " << via.id 
                 << ", Type: " << DieType2Str(via.type) << ", X: " << via.x << ", Y: " << via.y << ")\n";

            // **設定方向**
            string direction;
            if (i == 0) {  // 最右邊的 ViaNode => Upper
                direction = "top-right";
            } else if (i == 1) {  // 第二個 ViaNode => Lower
                direction = "bottom-right";
            } else {  // 其餘所有的 ViaNode
                direction = "top-right";
            }

            cout << "Routing Direction: " << direction << "\n\n";

            // **找最近的 EdgeNode**
            EdgeNode* current_edge = findNearestEdgeNode(via, RDL.edge_nodes, direction);
            if (!current_edge) {
                cout << "No suitable EdgeNode found.\n";
                continue;
            }
            
            // ** 取得 head 和 tail nodes**
            string head_node = current_edge->net_sequence.front();  // head node
            string tail_node = current_edge->net_sequence.back();   // tail node

            cout << "EdgeNode ID: " << current_edge->id << endl;
            cout << "Original Net Sequence: ";
            for (const auto& node : current_edge->net_sequence) {
                cout << node << " => ";
            }
            cout << "END" << endl;

            list<string> new_sequence;
            bool inserted = false;

            bool Upper = false; // 用於看 ground 要加在 signal 上還是下
            // ** 建立 `access`、`net` 和 `ground`**
            string access = "access";  // 插入 access
            string net = "Net_" + to_string(via.id);  // Net_viaX
            string ground = "ground";  // 插入 ground

            for (auto it = current_edge->net_sequence.begin(); it != current_edge->net_sequence.end(); ++it) {
                // **如果 `via` 是 `tail_node`，先檢查並插入**
                if (*it == tail_node && tail_node.substr(0, 3)=="SIG" && via.id == stoi(tail_node.substr(4))) { // => upper
                    auto prev_it = prev(it);  // 取得 tail_node 前一個元素
                    if (*prev_it != ground) { 
                        new_sequence.push_back(ground);  // **只在 tail_node 前面插入 ground (若沒有的話)**
                    }
                    new_sequence.push_back(net);
                    new_sequence.push_back(access);
                    inserted = true;
                    Upper = true;
                }
            
                // **加入原本的節點**
                new_sequence.push_back(*it);
            
                // **如果 `via` 是 `head_node`，則在 `head_node` 後面插入 `access`、`Net_via`、`ground`**
                if (*it == head_node && head_node.substr(0, 3)=="SIG" && via.id == stoi(head_node.substr(4))) { // => lower
                    new_sequence.push_back(access);
                    new_sequence.push_back(net);
                    auto next_it = next(it);  // 取得 head_node 後一個元素
                    if (*next_it != ground) { 
                        new_sequence.push_back(ground);  // **只在 head_node 後面插入 ground (若沒有的話)**
                    }
                    inserted = true;
                }
            }

            // **更新 `current_edge->net_sequence`**
            // -2 是要減去 head 和 tail node
            int newSize = new_sequence.size() - 2;
            // **檢查是否有 "access"，若有則再減 1** => access 不算一種 net 
            int access_count = count(new_sequence.begin(), new_sequence.end(), "access");
            newSize -= access_count;
            
            if(newSize <= current_edge->capacity){
                current_edge->net_sequence = new_sequence;
                for (auto& edge : RDL.edge_nodes) {
                    if (edge.id == current_edge->id) {
                        edge.net_sequence = new_sequence;  // **確保更新的是 RDL.edge_nodes 內的正確 EdgeNode**
                        break;
                    }
                }
            }else{
                cout<<"Capacity is not enough"<<endl;
                return;
            }
            

            // **輸出更新後的 Net Sequence**
            cout << "Updated Net Sequence: ";
            for (const auto& node : current_edge->net_sequence) {
                cout << node << " => ";
            }
            cout << "END" << endl;

            bool first_iteration = true;
            while (true) {
                cout << "********************\n";

                // **檢查 EdgeNode 連接的 ViaNodes**
                ViaNode via1 = current_edge->vias[0];
                ViaNode via2 = current_edge->vias[1];

                // **第一次進入迴圈時，跳過檢查，避免誤判**
                if (!first_iteration) {
                    if ((via1.type == via.type && via1.id == via.id) || (via2.type == via.type && via2.id == via.id)) {
                        cout << "Matching ViaNode found. Routing complete.\n";

                        string viaName = via.type + "_" + to_string(via.id);
                        
                        // **找到 net_sequence 的 front 和 back**
                        auto front_it = current_edge->net_sequence.begin();
                        auto next_front = next(front_it);
                        auto back_it = prev(current_edge->net_sequence.end()); // 指向最後一個元素
                        auto prev_back = prev(back_it);

                        bool modified = false;

                        // ** 檢查 `front` 和 `back`**
                        if (*front_it == viaName) {  // 檢查 `Net_viaX` 是否在 front
                            // cout<<"front"<<endl;
                            current_edge->net_sequence.insert(next(front_it), "access");  // **在 front_it 後面插入 access**
                            modified = true;
                            
                        }else if (*back_it == viaName) {  // 檢查 `Net_viaX` 是否在 back
                            // cout<<"back"<<endl;
                            if(*(prev_back) == "ground"){
                                current_edge->net_sequence.erase(prev_back);
                            }
                            current_edge->net_sequence.insert(back_it, "access");  // **在 back_it 前面插入 access**
                            modified = true;
                        }

                        // ** 更新 EdgeNode**
                        if (modified) {
                            for (auto& edge : RDL.edge_nodes) {
                                if (edge.id == current_edge->id) {
                                    edge.net_sequence = current_edge->net_sequence;  // **確保更新的是 RDL.edge_nodes 內的正確 EdgeNode**
                                    break;
                                }
                            }

                            // ** 輸出更新後的 Net Sequence**
                            cout << "Updated Net Sequence: ";
                            for (const auto& node : current_edge->net_sequence) {
                                cout << node << " => ";
                            }
                            cout << "END" << endl << endl;
                        }



                        break;
                    }
                }

                first_iteration = false;

                // ** 選擇角度最大的 CrossTileEdge**
                CrossTileEdge* next_edge = findLargestAngleCrossTileEdge(current_edge, RDL.cross_tile_edges, RDL);
                if (!next_edge) {
                    cout << "No valid CrossTileEdge found. Stopping.\n";
                    break;
                }

                // **更新當前 EdgeNode**
                current_edge = (next_edge->edges[0].id == current_edge->id) ? &next_edge->edges[1] : &next_edge->edges[0];

                for (auto& edge : RDL.edge_nodes) {
                    if (edge.id == current_edge->id) {
                        current_edge = &edge;  // **確保更新的是 RDL.edge_nodes 內的正確 EdgeNode**
                        break;
                    }
                }

                // **取得原 tail node**
                string tail_node = current_edge->net_sequence.back();  // 原本的 tail

                cout << "EdgeNode ID: " << current_edge->id << endl;
                cout << "Original Net Sequence: ";
                for (const auto& node : current_edge->net_sequence) {
                    cout << node << " => ";
                }
                cout << "END" << endl;

                // ** 建立 `net` 和 `ground`**
                string net = "Net_" + to_string(via.id);  // Net_viaX
                string ground = "ground";  // Ground 繞過 Net_via

                list<string> new_sequence;
                bool inserted = false;  // 確保只插入一次

                for (auto it = current_edge->net_sequence.begin(); it != current_edge->net_sequence.end(); ++it) {
                    // **如果當前節點是 tail_node，則在前面插入 `ground`、`Net_xx`、`ground`**
                    if (*it == tail_node && !inserted) {
                        // **確保左邊沒有 ground 才插入**
                        if (Upper && !new_sequence.empty() && new_sequence.back() != ground) {
                            new_sequence.push_back(ground);
                        }

                        new_sequence.push_back(net);

                        // **確保 `Net_xx` 和 `tail_node` 之間有 ground**
                        if ((row==rightmost_per_row.size()-1 && num_nodes==1) || (!Upper && *it != ground)) {
                            new_sequence.push_back(ground);
                        }

                        inserted = true;
                    }

                    // **加入原本的節點**
                    new_sequence.push_back(*it);
                }

                // **更新 `current_edge->net_sequence`**
                // -2 是要減去 head 和 tail node
                int newSize = new_sequence.size() - 2;
                // **檢查是否有 "access"，若有則再減 1** => access 不算一種 net 
                int access_count = count(new_sequence.begin(), new_sequence.end(), "access");
                newSize -= access_count;

                if(newSize <= current_edge->capacity){
                    current_edge->net_sequence = new_sequence;
                    for (auto& edge : RDL.edge_nodes) {
                        if (edge.id == current_edge->id) {
                            edge.net_sequence = new_sequence;  // **確保更新的是 RDL.edge_nodes 內的正確 EdgeNode**
                            break;
                        }
                    }
                }else{
                    cout<<"Capacity is not enough"<<endl;
                    return;
                }

                // **輸出更新後的 Net Sequence**
                cout << "Updated Net Sequence: ";
                for (const auto& node : current_edge->net_sequence) {
                    cout << node << " => ";
                }
                cout << "END" << endl << endl;

                
            }

            cout << "-------------------\n\n";
        }
        
    }
}

void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routingCoordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL){
    RoutingGraph RDL1, RDL2 ; 
    
    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> die1LeftmostBumps, die2LeftmostBumps ; // = FindLeftmostInEachRow(bumps);
    vector<Bump> die1DummyBumps, die2DummyBumps ; // = ProcessLeftAndRight(bumps, leftmostBumps, routingCoordinate, Hor_SPACING_X);
    vector<double> die1Coordinate, die2Coordinate ; 
    unordered_map<string, unordered_set<string>> net_group_map;

    ofstream outputViaFile, outputTriangulationFile, outputNetlistFile ; 
    //未來需要調整的部分!!!
    die1Coordinate = routingCoordinate ; die1Coordinate[2] = routingCoordinate[2]/2 - Hor_SPACING_X ; 
    die2Coordinate = routingCoordinate ; die2Coordinate[0] = routingCoordinate[2]/2 + Hor_SPACING_X ; 

    copy_if(bumps.begin(), bumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(bumps.begin(), bumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;

    die1LeftmostBumps = FindLeftmostInEachRow(die1Bumps) ;
    die2LeftmostBumps = FindLeftmostInEachRow(die2Bumps) ;

    die1DummyBumps = ProcessLeftAndRight(die1Bumps, die1LeftmostBumps, die1Coordinate, Hor_SPACING_X) ;
    die2DummyBumps = ProcessLeftAndRight(die2Bumps, die2LeftmostBumps, die2Coordinate, Hor_SPACING_X) ;

    // for(auto bump : die1Bumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    // for(auto bump : die1DummyBumps) RDL1.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    // for(auto bump : die2Bumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    // for(auto bump : die2DummyBumps) RDL2.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;

    // // 使用 CDT 進行三角化
    // RDL1.edge_nodes = Triangulation(RDL1.via_nodes, designRule);
    // RDL2.edge_nodes = Triangulation(RDL2.via_nodes, designRule);
    
    // FindAdjacentViaNodes(RDL1) ; //!!!!!
    // FindAdjacentViaNodes(RDL2) ; //!!!!!

    // AddAccessViaEdges(RDL1) ;
    // AddAccessViaEdges(RDL2) ;

    // AddCrossTileEdges(RDL1, designRule) ;
    // AddCrossTileEdges(RDL2, designRule) ;


    // //要改寫繞線順序!!!
    // auto [rightmost_per_row, remaining_per_row] = findUpperLowerOtherPerRow(RDL1.via_nodes, net_group_map);

        
    // routeUpperLower(rightmost_per_row, RDL1);  //RDL1 RDL2 分隔導致無法熱繞線

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


bool isValidTile(const EdgeNode &A, const EdgeNode &B, const EdgeNode &C) {
    // 確保 A, B, C 是不同的 EdgeNode
    if (A.id == B.id || A.id == C.id || B.id == C.id) return false;

    // 找出這三個 EdgeNodes 共同連接的 ViaNode
    vector<tuple<string, DieType, int>> via_list;

    // **將 A, B, C 連結的 ViaNode 加入 via_list**
    for (const auto &via : A.vias) via_list.push_back(tuple<string, DieType, int>{via.dieName, via.type, via.id});
    for (const auto &via : B.vias) via_list.push_back(tuple<string, DieType, int>{via.dieName, via.type, via.id});
    for (const auto &via : C.vias) via_list.push_back(tuple<string, DieType, int>{via.dieName, via.type, via.id});

    // **對 via_list 進行排序，確保相同的 ViaNode 會相鄰**
    sort(via_list.begin(), via_list.end());

    // **使用 unique 去除重複元素**
    via_list.erase(unique(via_list.begin(), via_list.end()), via_list.end());

    // **確保這三個 EdgeNodes 共享 3 個獨特的 ViaNodes**
    if (via_list.size() != 3) return false;

    return true;
}


tuple<ViaNode, ViaNode, ViaNode> findSharedViaNodes(const EdgeNode &A, const EdgeNode &B, const EdgeNode &C) {
    map<tuple<string, DieType, int>, ViaNode> via_map; // 用來確保 ViaNode 唯一性
    map<tuple<string, DieType, int>, int> via_count; // 記錄 ViaNode 出現次數

    // **計算 ViaNode 出現次數**
    for (const auto &via : A.vias) {
        via_map[{via.dieName, via.type, via.id}] = via;
        via_count[{via.dieName, via.type, via.id}]++;
    }
    for (const auto &via : B.vias) {
        via_map[{via.dieName, via.type, via.id}] = via;
        via_count[{via.dieName, via.type, via.id}]++;
    }
    for (const auto &via : C.vias) {
        via_map[{via.dieName, via.type, via.id}] = via;
        via_count[{via.dieName, via.type, via.id}]++;
    }

    // **篩選出所有共享的 ViaNodes**
    vector<ViaNode> shared_vias;
    for (const auto &[key, count] : via_count) {
        if (count == 2) { // 只加入出現兩次的 ViaNode
            shared_vias.push_back(via_map[key]);
        }
    }

    // **確保有 3 個共享的 ViaNode**
    if (shared_vias.size() != 3) {
        return {ViaNode(), ViaNode(), ViaNode()}; // 返回無效的 ViaNode
    }

    // **確保 shared_vias 對應到正確的 EdgeNode**
    ViaNode via_AB, via_BC, via_CA;
    for (const auto &via : shared_vias) {
        bool in_A = false, in_B = false, in_C = false;
        for (const auto &v : A.vias) {
            if (v.dieName == via.dieName && v.type == via.type && v.id == via.id) in_A = true;
        }
        for (const auto &v : B.vias) {
            if (v.dieName == via.dieName && v.type == via.type && v.id == via.id) in_B = true;
        }
        for (const auto &v : C.vias) {
            if (v.dieName == via.dieName && v.type == via.type && v.id == via.id) in_C = true;
        }

        if (in_A && in_B) via_AB = via;
        else if (in_B && in_C) via_BC = via;
        else if (in_C && in_A) via_CA = via;
    }

    return {via_AB, via_BC, via_CA};
}

double calculateAngle(const ViaNode &V, const EdgeNode &E1, const EdgeNode &E2) {
    // 向量 V → E1
    double vecA_x = E1.x - V.x;
    double vecA_y = E1.y - V.y;

    // 向量 V → E2
    double vecB_x = E2.x - V.x;
    double vecB_y = E2.y - V.y;

    // 計算內積
    double dot_product = (vecA_x * vecB_x) + (vecA_y * vecB_y);

    // 計算向量長度
    double lenA = sqrt(vecA_x * vecA_x + vecA_y * vecA_y);
    double lenB = sqrt(vecB_x * vecB_x + vecB_y * vecB_y);

    // **計算角度**
    double angle = acos(dot_product / (lenA * lenB)) * 180.0 / M_PI;  

    return angle;  // 回傳角度 (degree)
}

void AddCrossTileEdges(RoutingGraph &RDL, DesignRule designRule) {
    int cross_tile_id = static_cast<int>(RDL.cross_tile_edges.size());  // 記錄初始 ID
    unordered_set<string> added_tiles;  // 避免重複加入 Tile
    unordered_set<string> added_edges;  // 避免重複加入個別 Cross-Tile Edge
    
    // **遍歷所有 EdgeNode A**
    for (auto &edge_A : RDL.edge_nodes) {
        vector<EdgeNode*> potential_edges;
        
        // **找到所有透過 `AccessViaEdge` 連接到 EdgeNode A 的 ViaNodes**
        for (auto &access_via : RDL.access_via_edges) {
            if (access_via.edge.id == edge_A.id) {  
                ViaNode via = access_via.via;

                // **找到與此 ViaNode 連結的其他 EdgeNodes**
                for (auto &other_access_via : RDL.access_via_edges) {
                    if (other_access_via.edge.id != edge_A.id && 
                        other_access_via.via.dieName == via.dieName && 
                        other_access_via.via.type == via.type && 
                        other_access_via.via.id == via.id) 
                    {
                        // if(coutFlag<10){
                        //     ++coutFlag ; 
                        //     Net A("dummy", {tuple<double,double,double,double>{other_access_via.edge.x, other_access_via.edge.y, via.x, via.y}}) ;
                        //     debugNets.push_back(A) ;
                        // }
                        potential_edges.push_back(&other_access_via.edge);
                    }
                }
            }
        }

        // **檢查這些 EdgeNodes 是否能形成 Tile**
        for (size_t i = 0; i < potential_edges.size(); i++) {
            for (size_t j = i + 1; j < potential_edges.size(); j++) {
                EdgeNode* edge_B = potential_edges[i];
                EdgeNode* edge_C = potential_edges[j];

                // **檢查 Edge Node A, B, C 是否能構成 Tile**
                if (isValidTile(edge_A, *edge_B, *edge_C)) {
                    // **確保不重複加入該 Tile**
                    string tile_key = to_string(edge_A.id) + "-" + to_string(edge_B->id) + "-" + to_string(edge_C->id);
                    if (added_tiles.find(tile_key) == added_tiles.end()) {
                        // **找出共享的 ViaNodes**
                        auto [via_AB, via_BC, via_CA] = findSharedViaNodes(edge_A, *edge_B, *edge_C);
                        double angle_AB = calculateAngle(via_AB, edge_A, *edge_B);
                        double angle_BC = calculateAngle(via_BC, *edge_B, *edge_C);
                        double angle_CA = calculateAngle(via_CA, *edge_C, edge_A);

                        // **計算長度**
                        double length_AB = sqrt(pow(edge_A.x - edge_B->x, 2) + pow(edge_A.y - edge_B->y, 2));
                        double length_BC = sqrt(pow(edge_B->x - edge_C->x, 2) + pow(edge_B->y - edge_C->y, 2));
                        double length_CA = sqrt(pow(edge_C->x - edge_A.x, 2) + pow(edge_C->y - edge_A.y, 2));


                        // ** 逐個加入 Cross-Tile Edge，確保不重複**
                        string edge_AB_key = to_string(edge_A.id) + "-" + to_string(edge_B->id);
                        string edge_BC_key = to_string(edge_B->id) + "-" + to_string(edge_C->id);
                        string edge_CA_key = to_string(edge_C->id) + "-" + to_string(edge_A.id);

                        if (added_edges.find(edge_AB_key) == added_edges.end()) {
                            RDL.cross_tile_edges.emplace_back(cross_tile_id++, vector<EdgeNode>{edge_A, *edge_B}, angle_AB);
                            added_edges.insert(edge_AB_key);
                        }
                        if (added_edges.find(edge_BC_key) == added_edges.end()) {
                            RDL.cross_tile_edges.emplace_back(cross_tile_id++, vector<EdgeNode>{*edge_B, *edge_C}, angle_BC);
                            added_edges.insert(edge_BC_key);
                        }
                        if (added_edges.find(edge_CA_key) == added_edges.end()) {
                            RDL.cross_tile_edges.emplace_back(cross_tile_id++, vector<EdgeNode>{*edge_C, edge_A}, angle_CA);
                            added_edges.insert(edge_CA_key);
                        }

                        added_tiles.insert(tile_key);  // 避免重複
                    }
                }
            }
        }
    }
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
        
            
        string src_type = DieType2Str(point_to_via[segment.source()].type) + "_" + to_string(point_to_via[segment.source()].id);
        string tgt_type = DieType2Str(point_to_via[segment.target()].type) + "_" + to_string(point_to_via[segment.target()].id);

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
