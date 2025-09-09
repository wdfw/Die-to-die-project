#pragma once

#include <filesystem>
#include <vector>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DesignRule.hpp"
#include "Bump.hpp"
#include "Utils.hpp"
#include "RoutingGraph.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

extern vector<Bump> debugBumps ; 

using namespace std ; 

int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routingCoordinate, 
                const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes) ;

void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routingCoordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL) ; 

vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& routingCoordinate, double Hor_SPACING_X) ;
vector<Bump> FindLeftmostInEachRow(const vector<Bump>& bumps) ;
void FindAdjacentViaNodes(RoutingGraph &RDL) ; 

vector<EdgeNode> Triangulation(const vector<ViaNode>& via_nodes, const DesignRule& designRule) ;
vector<EdgeNode*> FindAdjacentEdgeNodes( ViaNode &via,  vector<EdgeNode> &edge_nodes) ;

vector<EdgeNode> CreateEdgeNodes(const vector<EdgeNode>& trigulationNodes) ; 
void AddAccessViaEdges(RoutingGraph &RDL) ;

class D2DRouting{
private:
    vector<EdgeNode> Triangulation(const vector<ViaNode>& via_nodes, const DesignRule& designRule) ;
    vector<Bump> FindLeftmostInEachRow(const vector<Bump>& bumps) ;
    vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& routingCoordinate, double Hor_SPACING_X) ;

    void AddAccessViaEdges(RoutingGraph &RDL) ;
    void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routingCoordinate, const DesignRule& designRule, 
                               const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL) ; \
public:
    int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routingCoordinate, 
                const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes) ;

} ;