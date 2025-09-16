#include "RoutingGraph2.hpp"

//------------------------------------------ ViaNode2 Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const ViaNode2& node) {
    // os << *(static_cast<const Bump*>(node)) << ":\n" ; 
    for(auto& tileNode : node.tileNodes){
        // os << static_cast<const Bump&>(*tileNode) << ", " ; 
    }
    return os ; 
}

//------------------------------------------ ViaNode2 Method End ------------------------------------------ 

//------------------------------------------ TileNode2 Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const TileNode2& node) {
    // os << static_cast<const Bump&>(*node) << ":\n" ; 
    // os << "Tiles: " ;
    // for(auto& tileNode : node.tileNodes) os << static_cast<const Bump&>(*tileNode) << " " <<  << ", " ; 

    // os << "\nVias: " ;
    // for(auto& viaNode : node.viaNodes) os << static_cast<const Bump&>(*viaNode) << ", " ; 
        
    return os ; 
}

//------------------------------------------ TileNode2 Method End ------------------------------------------ 
