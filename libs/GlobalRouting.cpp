#include "GlobalRouting.hpp"

const double EPSILON_Y = 1e-6 ;
const double EPSILON_X = 1e-3;

int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routing_area_coordinate, 
                const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes){

    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> marginalBumps ;
    
    double min_y = numeric_limits<double>::max() ;
    double Hor_SPACING_X ; 
    
    string dirPath ; 
    
    Timer timer;

    for(const auto& bump : bumps){
        if(bump.type == SIGNAL){
            if(bump.y<min_y){
                min_y = bump.y ;
                marginalBumps.clear() ; 
            }
            if(fabs(bump.y - min_y) < EPSILON_Y){
               marginalBumps.push_back(bump) ;
            }
        }
    }
    sort(marginalBumps.begin(), marginalBumps.end(), [](const Bump& a, const Bump& b) {return a.x < b.x;}) ;

    Hor_SPACING_X = marginalBumps[1].x - marginalBumps[0].x ; 

    for(int layer = 1, leftBumpCount = bumps.size(); leftBumpCount; ++layer){
        dirPath = output_path + "RDL" + to_string(layer) + "/";

        if (!filesystem::exists(dirPath)) filesystem::create_directories(dirPath);

        
        timer.SetClock() ;
        
        ConstructRoutingGraph(bumps, routing_area_coordinate, designRule, dirPath, Hor_SPACING_X, layer, allRDL);

        globalRouteTimes.push_back(timer.GetDurationMilliseconds()); // **記錄每層的 Global Route 時間**
        leftBumpCount = 0 ;
        // for(const auto& bump : die1Bumps) if(bump.type==SIGNAL) ++leftBumpCount ;

    }

    return globalRouteTimes.size() ; 
}

void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routing_area_coordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL){
    RoutingGraph RDL; 
    
    vector<Bump> die1Bumps, die2Bumps ; 
    vector<Bump> die1LeftmostBumps, die2LeftmostBumps ; // = FindLeftmostInEachRow(bumps);
    vector<Bump> die1DummyBumps, die2DummyBumps ; // = ProcessLeftAndRight(bumps, leftmostBumps, routing_area_coordinate, Hor_SPACING_X);
    vector<double> die1Coordinate, die2Coordinate ; 

    //
    die1Coordinate = routing_area_coordinate ; die1Coordinate[2] = routing_area_coordinate[2]/2 - Hor_SPACING_X ; 
    die2Coordinate = routing_area_coordinate ; die2Coordinate[0] = routing_area_coordinate[2]/2 + Hor_SPACING_X ; 

    copy_if(bumps.begin(), bumps.end(), back_inserter(die1Bumps), [](const Bump& bump){return bump.name=="DIE1"; }) ;
    copy_if(bumps.begin(), bumps.end(), back_inserter(die2Bumps), [](const Bump& bump){return bump.name=="DIE2"; }) ;

    die1LeftmostBumps = FindLeftmostInEachRow(die1Bumps) ;
    die2LeftmostBumps = FindLeftmostInEachRow(die2Bumps) ;

    die1DummyBumps = ProcessLeftAndRight(die1Bumps, die1LeftmostBumps, die1Coordinate, Hor_SPACING_X) ;
    die2DummyBumps = ProcessLeftAndRight(die2Bumps, die2LeftmostBumps, die2Coordinate, Hor_SPACING_X) ;

    for(auto& it : die1Bumps) cout << it << "\n" ;
    for(auto& it : die2Bumps) cout << it << "\n" ;
    for(auto& it : die1DummyBumps) cout << it << "\n" ;
    for(auto& it : die2DummyBumps) cout << it << "\n" ;

    // for(auto bump : bumps){
    //     RDL.via_nodes.emplace_back(bump.name, bump.type, bump.id, bump.x, bump.y) ;
    // }
        
    // for(auto dummy : dummys){
    //     RDL.via_nodes.push_back("dum", DUMMY, dummy.id, dummy.x, dummy.y);
    // }

    // 使用 CDT 進行三角化
    // RDL.edge_nodes = triangulation(RDL.via_nodes, output_path, designRule);
    

}

vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& coordinate, double Hor_SPACING_X){
    double min_x = coordinate[0], max_x = coordinate[2] ;

    Bump newPoint ;
    vector<Bump> dummies;

    double current_x ; 
    
    for (const auto& point : leftmostBumps) {
        double current_x = point.x;

        // 向左移動
        while (current_x - Hor_SPACING_X >= min_x) {
            current_x -= Hor_SPACING_X;
            newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
                dummies.push_back(newPoint);
            }
        }

        // 向右移動
        current_x = point.x;
        while (current_x + Hor_SPACING_X <= max_x) {
            current_x += Hor_SPACING_X;
            newPoint = Bump("Dummy", DUMMY, dummies.size(), current_x, point.y);

            if(find_if(bumps.begin(), bumps.end(), [&newPoint](const Bump& bump){return fabs(bump.x - newPoint.x) < EPSILON_X && fabs(bump.y - newPoint.y) < EPSILON_Y;})==bumps.end()){
                dummies.push_back(newPoint);
            }
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

