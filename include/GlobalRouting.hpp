#pragma once

#include <filesystem>
#include <vector>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DesignRule.hpp"
#include "Bump.hpp"
#include "Utils.hpp"
#include "RoutingGraph.hpp"
#include "RoutingGraph2.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

extern vector<Bump> debugBumps ; 


using namespace std ; 

int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routingCoordinate, 
                const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes) ;

vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& routingCoordinate, double Hor_SPACING_X) ;
vector<Bump> FindLeftmostInEachRow(const vector<Bump>& bumps) ;
void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routingCoordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL) ; 

vector<EdgeNode> Triangulation(const vector<ViaNode>& via_nodes, const DesignRule& designRule) ;

void FindAdjacentViaNodes(RoutingGraph &RDL) ; 

vector<EdgeNode*> FindAdjacentEdgeNodes( ViaNode &via,  vector<EdgeNode> &edge_nodes) ;
void AddAccessViaEdges(RoutingGraph &RDL) ;

vector<EdgeNode> CreateEdgeNodes(const vector<EdgeNode>& trigulationNodes) ; 
void AddCrossTileEdges(RoutingGraph &RDL, DesignRule designRule) ;
void AddCrossTileEdges2(RoutingGraph &RDL, DesignRule designRule) ;

pair<vector<vector<ViaNode>>, vector<vector<ViaNode>>> findUpperLowerOtherPerRow(const vector<ViaNode>& vias, unordered_map<string, unordered_set<string>>& net_group_map) ;
void routeUpperLower(vector<vector<ViaNode>>& rightmost_per_row, RoutingGraph& RDL) ;
EdgeNode* findNearestEdgeNode(const ViaNode& via, vector<EdgeNode>& edge_nodes, const string& direction) ;
CrossTileEdge* findLargestAngleCrossTileEdge(EdgeNode* current_edge, vector<CrossTileEdge>& cross_tile_edges, RoutingGraph& RDL) ;