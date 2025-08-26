#include "Parser.hpp"

void ParseBump(const string &inputPath, vector<Bump> &bumps, vector<double> &coordinates){ 
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseBump] Failed to open bump file: " + inputPath);

    string line, dieName, type ;
    int id ;
    double x, y, v1, v2 ; ;

    bumps.clear() ; coordinates.clear() ; 

    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;
        istringstream iss(line);
        
        if(lineNumber<=2 && (iss >> v1 >> v2)){
            coordinates.push_back(v1) ; coordinates.push_back(v2) ;
        }else{
            if (!(iss >> dieName >> type >> id >> x >> y )) throw runtime_error("[ParseBump] Error parsing in line #" + to_string(lineNumber)) ;
            bumps.emplace_back(dieName, Str2DieType(type), id, x, y);
        }
    }
}

void ParseDesignRule(const string &inputPath, DesignRule &designRule){
    
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseDesignRule] Failed to open design rule file: " + inputPath);

    string line;
    
    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;

        istringstream iss(line);
        string key;
        if (!(iss >> key)) throw runtime_error("[ParseDesignRule] Error parsing in line #" + to_string(lineNumber)) ;
        if(key == "via_opening_diameter"){
            iss >> designRule.viaOpeningDiameter ;
        }else if(key == "via_pad_diameter"){
            iss >> designRule.viaPadDiameter ;
        }else if(key == "minimum_via_pad_spacing"){
            iss >> designRule.minimumViaPadSpacing ;
        }else if(key == "minimum_via_spacing"){
            iss >> designRule.minimumViaSpacing ; 
        }else if(key == "minimum_line_width"){
            iss >> designRule.minimumLineWidth ;
        }else if(key == "minimum_line_spacing"){
            iss >> designRule.minimumLineSpacing ;
        }else if(key == "minimum_line_via_spacing"){
            iss >> designRule.minimumLineViaSpacing ;
        }else if(key == "minimum_teardropDist"){
            iss >> designRule.minimumTeardropDist ;
        }else{
            cerr << "Unknown Rule:" << key << endl;
        }
    }
    if(designRule.viaOpeningDiameter>designRule.viaPadDiameter){
        throw runtime_error("[ParseDesignRule] Invalid design rule : via opening diameter larger than via pad diameter");
    }
    file.close();
}

void ParseOffsetBump(const string &inputPath, vector<Bump> &bumps, vector<Bump> &matchedBumps){
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseOffsetBump] Failed to open offset bump file: " + inputPath);

    string line, dieName, type ;
    int id ;
    double x1, y1, x2, y2 ; ;

    bumps.clear() ; matchedBumps.clear() ;
    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;
        istringstream iss(line);
        
        if (!(iss >> dieName >> type >> id >> x1 >> y1 >> x2 >> y2 )) throw runtime_error("[ParseOffsetBump] Error parsing in line #" + to_string(lineNumber)) ;
        // if( matchedBumps.find(matches[{dieName, Str2DieType(type), id}])==matchedBumps.end() ) throw runtime_error("[ParseOffsetBump] Unmatch bum in line #" + to_string(lineNumber)) ;

        bumps.emplace_back(dieName, Str2DieType(type), id, x1, y1) ;
        matchedBumps.emplace_back(dieName, Str2DieType(type), id, x2, y2) ;
    }
}

void ParseNet(const string &inputPath, vector<Net> &nets){
    
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseNet] Failed to open netlist file: " + inputPath);

    string line, netName ;
    map<string, int> usedName ;

    double x1, y1, x2, y2 ; ;

    nets.clear() ; 
    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;
        istringstream iss(line);
        if (!(iss >> netName >> x1 >> y1 >> x2 >> y2 )) throw runtime_error("[ParseNet] Error parsing in line #" + to_string(lineNumber)) ;
        if(usedName.find(netName)==usedName.end()){
            usedName[netName] = nets.size() ;
            nets.emplace_back(Net(netName)) ;
        }
        nets[usedName[netName]].emplace_back(x1, y1, x2, y2) ;
    }
}

void ParseTeardrop(const string &inputPath, vector<tuple<Bump, double, double, double, double>> &teardrops){
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseTeardrop] Failed to open teardrop file: " + inputPath);


    string line, dieName, type ;
    int id ;
    double x1, y1, x2, y2 ; ;

    teardrops.clear() ; 
    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;

        istringstream iss(line);
        if (!(iss >> dieName >> type >> id >> x1 >> y1 >> x2 >> y2 )) throw runtime_error("[ParseTeardrop] Error parsing in line #" + to_string(lineNumber)) ;
        teardrops.push_back({{dieName, Str2DieType(type), id, x1, y1}, x1, y1, x2, y2}) ; //x1, y1是bump位置, x2, y2是目標線段起點的位置
    }
}

void ParseTriangulation(const string &inputPath, vector<tuple<double, double, double, double>> &triangulationEdges){
    ifstream file(inputPath);
    if(!file.is_open()) throw runtime_error("[ParseTriangulation] Failed to open triangulation file: " + inputPath);

    string line, dieName, type ;
    double x1, y1, x2, y2 ; ;

    triangulationEdges.clear() ; 
    for(int lineNumber=1; getline(file, line); ++lineNumber){
        if((line = Strip(line)).empty()) continue ;

        istringstream iss(line);
        if (!(iss >> x1 >> y1 >> x2 >> y2 )) throw runtime_error("[ParseTriangulation] Error parsing in line #" + to_string(lineNumber)) ;
        triangulationEdges.push_back({x1, y1, x2, y2}) ; 
    }
}
