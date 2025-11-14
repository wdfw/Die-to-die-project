#pragma once
#include <vector>
#include <map>
#include <set>
#include <memory>

#include "Bump.hpp"

class ViaNode2 ; 
class TileNode2 ; 
class RoutingGraph2 ; 
struct TileToTileEdge ;

ostream& operator<<(ostream& os, const ViaNode2& node)  ;
ostream& operator<<(ostream& os, const TileNode2& node) ;


enum EdgeType {
    Undefined,
    BaseEdge, 
    LegEdge
} ; 

struct EdgeNode {
    shared_ptr<ViaNode2> viaNode1, viaNode2 ; 
    EdgeType type ; 
    vector<pair<double, double>> positions ; 
    EdgeNode(shared_ptr<ViaNode2> viaNode1, shared_ptr<ViaNode2> viaNode2, EdgeType type, const vector<pair<double, double>> positions) :
        viaNode1(viaNode1), viaNode2(viaNode2), type(type), positions(positions) {}
} ; 

struct TileToTileEdge : public shared_ptr<TileNode2> {
    shared_ptr<int> capacity = make_shared<int>(0) ;
    shared_ptr<int> currentCapacity = make_shared<int>(0) ;
    shared_ptr<ViaNode2> crossedViaNode1, crossedViaNode2  ; 
    TileToTileEdge() = default ; 
    TileToTileEdge(const shared_ptr<TileNode2>& node) : shared_ptr<TileNode2>(node) {} ;
} ;

struct GlobalNet : public vector<TileToTileEdge>{
    string name;
    DieType type;
    int id;

    shared_ptr<ViaNode2> startViaNode ; 
    shared_ptr<ViaNode2> endViaNode ; 

    GlobalNet(const string name="", const DieType& type=DUMMY, int id=-1, shared_ptr<ViaNode2> startViaNode=nullptr, shared_ptr<ViaNode2> endViaNode=nullptr)
                : name(name), type(type), id(id), startViaNode(startViaNode), endViaNode(endViaNode) {}

};

struct DetailedNet : public vector<pair<shared_ptr<EdgeNode>, int>>{
    string name;
    DieType type;
    int id;

    shared_ptr<ViaNode2> startViaNode ; 
    shared_ptr<ViaNode2> endViaNode ; 

    DetailedNet(const string name="", const DieType& type=DUMMY, int id=-1, shared_ptr<ViaNode2> startViaNode=nullptr, shared_ptr<ViaNode2> endViaNode=nullptr)
                : name(name), type(type), id(id), startViaNode(startViaNode), endViaNode(endViaNode) {}

};

class ViaNode2 : public Bump {
public:
    vector<shared_ptr<ViaNode2>> viaNodes ; 
    vector<TileToTileEdge> tileNodes ; 
    using Bump::Bump ; 
    ViaNode2(const Bump& bump) : Bump(bump) {} ;
};

class TileNode2 : public Bump {
public:
    vector<shared_ptr<ViaNode2>> viaNodes ; 
    vector<TileToTileEdge> tileNodes ; 
    using Bump::Bump ; 
    TileNode2(const Bump& bump) : Bump(bump) {} ;
};

// class Node2 : public Bump {
// public:
//     vector<shared_ptr<ViaNode2>> viaNodes ; 
//     vector<TileToTileEdge> tileNodes ; 

//     using Bump::Bump ; 
//     TileNode2(const Bump& bump) : Bump(bump) {} ;
// };
class RoutingGraph2 {
public:
    vector<shared_ptr<ViaNode2>> viaNodes ;             
    vector<shared_ptr<TileNode2>> tileNodes ;       
    vector<shared_ptr<EdgeNode>> edgeNodes ;
};


