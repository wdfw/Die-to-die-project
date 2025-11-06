#pragma once

#include <filesystem>
#include <vector>
#include <limits>

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
extern vector<Net> debugNets ; 
extern vector<tuple<string, double, double>> debugLabels ; 
using namespace std ; 

class Router {
protected:
    const double epsilonX = 1e-3 ;
    const double epsilonY = 1e-6 ;
    void Initialize() ;
    // void FindLeftMostInEachRow(const vector<Bump>& bumps, vector<Bump>& leftMostBumps) ;

    void FindHorizontalSpace(const vector<Bump>& bumps, double& minimumHorizontalSpace) ; 
    void PaddingBumps(const vector<Bump>& bumps, vector<Bump>& paddingBumps, const vector<double>& coordinate, double minimumHorizontalSpace) ; 
    
    void CreateViaNodes(RoutingGraph2& graph, const vector<Bump>& bumps) ;
    void Triangulation(RoutingGraph2& graph) ;
    void ConnectTileTileEdges(RoutingGraph2& graph) ;
    void SetCapacity(RoutingGraph2& graph, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps) ; 

    void CreateEdgeNodes(RoutingGraph2& graph) ;

    void CombineRDLs(RoutingGraph2& graph1, RoutingGraph2& graph2, RoutingGraph2& graph3, const vector<double>& coordinate1, const vector<double>& coordinate2) ; 

    void SelectRoutingBumps(const vector<Bump>& bumps, vector<Bump>& routingBumps, vector<Bump>& offsetBumps, int selectNum = numeric_limits<int>::max()) ; 
    void CreateViaBumps(const vector<Bump>& offsetBumps, vector<Bump>& viaBumps) ;
    void ConstructRoutingGraph(const vector<Bump>& routingBumps, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps, RoutingGraph2& graph, double minimumHorizontalSpace, const vector<double>& coordinate) ; 
    void GenerateGraphFile(const vector<Bump>& routingBumps, const vector<Bump>& offsetBumps, const vector<Bump>& viaBumps, const RoutingGraph2& graph, int layer,  const string& directoryPath) ;
    
    virtual double GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph, vector<GraphNet>& nets) ; 
    
    virtual double DetailRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph) ; 
public:
    DesignRule designRule ; 
    vector<Bump> bumps ; 
    vector<RoutingGraph2> routingGraphs ;
    vector<double> coordinate ;
    vector<clock_t> routingTimes ;

    Router(const DesignRule& designRule, const vector<Bump>& bumps, const vector<double>& coordinate) ;
    void Solve(const string& outputDirectories="") ; 
} ;

