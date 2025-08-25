#include "Bump.hpp"

//------------------------------------------ Bump Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const Bump& bump){
    os << bump.name << " " << DieType2Str(bump.type) << " " << bump.id << " " << bump.x << " " << bump.y ;
    return os ; 
}

//------------------------------------------ Bump Method End ------------------------------------------ 

//------------------------------------------ Net Method Begin ------------------------------------------ 

ostream& operator<<(ostream& os, const Net& bump){
    os << bump.name << ": " ;
    for(auto& it : bump) os << "(" << get<0>(it) << "," << get<1>(it) << "," << get<2>(it) << "," << get<3>(it) << ")" ;
    return os ;
}

//------------------------------------------ Net Method Begin ------------------------------------------ 

string DieType2Str(const DieType& type){
    switch (type){
        case DUMMY : return "Dummy" ;
        case SIGNAL :  return "Signal" ;
        case VDD :  return "Vdd" ;
        case VSS :  return "Vss" ;
    }
    return "UNKNOWN" ;
}

DieType Str2DieType(const string& str){
    if(str=="Dummy") return DUMMY ;
    else if(str=="Signal" || str=="SIG") return SIGNAL ;
    else if(str=="Vdd" || str=="VDD") return VDD ;
    else if(str=="Vss" || str=="VSS") return VSS ;
    throw runtime_error("Can't convert string \"" + str + "\" to DieType") ;
}
