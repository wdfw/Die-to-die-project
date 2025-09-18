#pragma once
#include <vector>
#include <map>
#include <set>
#include <memory>
#include "Bump.hpp"

class ViaNode2 ; 
class TileNode2 ; 
class RoutingGraph2 ; 

ostream& operator<<(ostream& os, const ViaNode2& node)  ;
ostream& operator<<(ostream& os, const TileNode2& node) ;

class ViaNode2 : public Bump {
public:
    vector<shared_ptr<ViaNode2>> viaNodes ; 
    vector<shared_ptr<TileNode2>> tileNodes ; 
    using Bump::Bump ; 
    ViaNode2(const Bump& bump) : Bump(bump) {} ;
};

class TileNode2 : public Bump {
public:
    struct TileToTileEdge : public shared_ptr<TileNode2> {
        shared_ptr<int> capacity = make_shared<int>(0) ;
        // shared_ptr<pair<shared_ptr<TileNode2>, shared_ptr<TileNode2>>> crossedEdge ;
        TileToTileEdge(const shared_ptr<TileNode2>& node) : shared_ptr<TileNode2>(node) {} ;
    } ;
    vector<shared_ptr<ViaNode2>> viaNodes ; 
    vector<TileToTileEdge> tileNodes ; 

    using Bump::Bump ; 
    TileNode2(const Bump& bump) : Bump(bump) {} ;
};

class RoutingGraph2 {
public:
    vector<shared_ptr<ViaNode2>> viaNodes;             
    vector<shared_ptr<TileNode2>> tileNodes;            
    map<shared_ptr<ViaNode2>, set<shared_ptr<TileNode2>>> viaToTileEdges ; 
    map<shared_ptr<TileNode2>, set<shared_ptr<TileNode2>>> tileToTileEdges ; 


};

