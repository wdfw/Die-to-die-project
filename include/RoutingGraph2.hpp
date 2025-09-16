#pragma once
#include <vector>
#include <map>
#include <set>
#include "Bump.hpp"

class ViaNode2 ; 
class TileNode2 ; 
class RoutingGraph2 ; 

ostream& operator<<(ostream& os, const ViaNode2& node)  ;
ostream& operator<<(ostream& os, const TileNode2& node) ;

class ViaNode2 : public Bump {
public:
    vector<TileNode2*> tileNodes ; 

    using Bump::Bump ; 
    ViaNode2(const Bump& bump) : Bump(bump) {} ;
};

class TileNode2 : public Bump {
    struct TileToTileEdge{
        TileNode2* node ; 
        int capicity ;
    } ;
public:
    vector<ViaNode2*> viaNodes ; 
    vector<TileToTileEdge> tileNodes ; 

    using Bump::Bump ; 
    TileNode2(const Bump& bump) : Bump(bump) {} ;
};


class RoutingGraph2 {
public:
    vector<ViaNode2> viaNodes;             
    vector<TileNode2> tileNodes;            
   
    map<ViaNode2*, set<TileNode2*>> viaToTileEdges ; 
    map<TileNode2*, set<TileNode2*>> tileToTileEdges ; 
};

