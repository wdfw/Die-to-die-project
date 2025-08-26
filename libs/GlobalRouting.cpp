#include "GlobalRouting.hpp"

int RoutingGuideGenerate(vector<Bump> vias, vector<double> routing_area_coordinate, DesignRule designRule, string output_path, vector<RoutingGraph>& allRDL, double &offset, vector<int>& GlobalRouteTime){
    double min_y = numeric_limits<double>::max();
    vector<Bump> filtered_vias;

    // if(offset == 0)
    //     EPSILON_Y = 1e-6;
    // else
    //     EPSILON_Y = offset * 10 + 5;

    // **找出 DIE1 中 y 值最小的 vias**
    for (const auto& via : vias) {
        if (via.dieName == "DIE1" && via.type == "SIG") {
            if (via.y < min_y) {
                min_y = via.y;
            }
        }
    }

    // **過濾出 y 值最小的 vias**
    for (const auto& via : vias) {
        if (via.dieName == "DIE1" && via.type == "SIG" && fabs(via.y - min_y) < 1e-6) {
            filtered_vias.push_back(via);
        }
    }

    // **按照 x 值排序**
    sort(filtered_vias.begin(), filtered_vias.end(), [](const Bump& a, const Bump& b) {
        return a.x < b.x;
    });
    
    double Hor_SPACING_X = filtered_vias[1].x - filtered_vias[0].x; // 找到兩個 bump 之間的距離, 用來抓 dummy point

    int layer = 1;
    string dirPath = output_path;
    GlobalRouteTime.reserve(100); // 預留空間，避免多次 reallocate
    while(true){
        int sigSize = 0;
        for (const auto& via : vias) {
            cout << "(" << via.type << ", " << via.id << ", " << via.dieName << ") " << endl;
            if(via.type == "SIG")
                sigSize++;
        }
        if(sigSize == 0){
            cout << "\nGlobal Route: all net has been routed\n\n";
            return layer - 1; // 因為 layer 多加一次
        }else{
            dirPath = output_path + "RDL" + to_string(layer) + "/";
            // 建立資料夾（如果不存在）
            if (!fs::exists(dirPath)) {
                fs::create_directories(dirPath);
            }

            cout<< "Routing layer " << layer << endl << endl;

            if(layer != 1){
                if(layer == 2){ // 從 RDL 2 開始，最左側的兩顆 VSS via 不會繼續往上拉線, 因為做 via offset 時會檔到
                    vector<Bump> vss_vias;
                    for (auto const& via : vias) {
                        if (via.type == "VSS") {
                            vss_vias.push_back(via);
                        }
                    }
                    sort(vss_vias.begin(), vss_vias.end(),
                        [](auto const& a, auto const& b) {
                            return a.x < b.x;
                    });

                    // 取出最小的兩個 id
                    int id1 = vss_vias[0].id;
                    int id2 = vss_vias[1].id;

                    for (auto &via : vias) {
                        if (via.type == "VSS" && (via.id == id1 || via.id == id2)) {
                            via.type = "dum";
                        }
                    }

                    routing_area_coordinate[0] -= 100;
                }   

                ofstream outFile(dirPath + "via_layer_" + to_string(layer));
                outFile << routing_area_coordinate[0] << " " << routing_area_coordinate[1] << "\n";
                outFile << routing_area_coordinate[2] << " " << routing_area_coordinate[3] << "\n";
                
                for (auto& via : vias) {
                    if(layer != 2)
                        via.x -= (designRule.minimum_via_spacing + designRule.via_opening_diameter); // 因為 via offset 會讓 via 的 x 座標偏移，所以這邊要先做 offset
                    if(via.type == "SIG" || via.type == "VSS")
                        outFile << via.dieName << " " << via.type << " " << via.id << " " << via.x << " " << via.y << endl;
                    if(via.type == "VDD") // 為了不要讓 power 出現在 RDL 2,3,4
                       via.type = "dum";
                }
                outFile.close();
            }

            Timer timer;
            // constructRoutingGraph(vias, routing_area_coordinate, designRule, dirPath, Hor_SPACING_X, layer, allRDL);
            
            GlobalRouteTime.push_back(timer.GetDurationMilliseconds()); // **記錄每層的 Global Route 時間**
            
            layer++;
        }
    }
    
}