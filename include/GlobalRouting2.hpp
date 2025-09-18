#pragma once

#include <filesystem>
#include <vector>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DesignRule.hpp"
#include "Bump.hpp"
#include "Utils.hpp"
#include "RoutingGraph2.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

extern vector<Bump> debugBumps ; 

using namespace std ; 

// int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routingCoordinate, 
//                 const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
//                 double offset, vector<clock_t>& globalRouteTimes) ;

// vector<EdgeNode*> FindAdjacentEdgeNodes( ViaNode &via,  vector<EdgeNode> &edge_nodes) ;
// void AddAccessViaEdges(RoutingGraph &RDL) ;

// vector<EdgeNode> CreateEdgeNodes(const vector<EdgeNode>& trigulationNodes) ; 
// void AddCrossTileEdges(RoutingGraph &RDL, DesignRule designRule) ;
// void AddCrossTileEdges2(RoutingGraph &RDL, DesignRule designRule) ;

class Router {
private:
    const double epsilonX = 1e-3 ;
    const double epsilonY = 1e-6 ;
    
    void FindLeftmostInEachRow(const vector<Bump>& bumps, vector<Bump>& leftMostBumps) ;
    void PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, const vector<double>& coordinate, double minimumHorizontalSpace) ; //生成 Padding Bump使其可以順利三角化
    
    void CreateViaNodes(RoutingGraph2& graph, const vector<Bump>& bumps) ;
    void Triangulation(RoutingGraph2& graph) ;
    void ConnectTileTileEdges(RoutingGraph2& graph) ;

    void CombineDie1Die2(RoutingGraph2& graph1, RoutingGraph2& graph2, const vector<double>& coordinate1, const vector<double>& coordinate2) ; 

    void ConstructRoutingGraph(const string& outputPath, int layer, double minimumHorizontalSpace) ; 
public:
    DesignRule designRule ; 
    vector<Bump> bumps ; 
    vector<RoutingGraph2> routingGraphs ;
    vector<double> coordinate ;
    vector<clock_t> routingTimes ;

    Router(const DesignRule& designRule, const vector<Bump>& bumps, const vector<double>& coordinate) ;

    void Initial() ;
    void GlobalRoute(const string& outputDirectories="") ; 
} ;