#include "RoutingGraph2.hpp"

//------------------------------------------ ViaNode2 Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const ViaNode2& node) {
    os << DieType2Str(node.type) << " " << node.id  ;
    // os << DieType2Str(node.type) << " " << node.id << ":\n" ;
    // os << "Tiles: " ;
    // for(auto& tileNode : node.tileNodes) os << DieType2Str(tileNode->type) << " " << tileNode->id << ", " ;

    // os << "\nVias: " ;
    // for(auto& viaNode : node.viaNodes) os << DieType2Str(viaNode->type) << " " << viaNode->id << ", " ;

    return os ; 
}

//------------------------------------------ ViaNode2 Method End ------------------------------------------ 

//------------------------------------------ TileNode2 Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const TileNode2& node) {
    os << DieType2Str(node.type) << " " << node.id  ;

    // os << DieType2Str(node.type) << " " << node.id << ":\n" ;

    // os << "Tiles: " ;
    // for(auto& tileNode : node.tileNodes) os << DieType2Str(tileNode->type) << " " << tileNode->id << ", " ;

    // os << "\nVias: " ;
    // for(auto& viaNode : node.viaNodes) os << DieType2Str(viaNode->type) << " " << viaNode->id << ", " ;
        
    return os ; 
}

//------------------------------------------ TileNode2 Method End ------------------------------------------ 
