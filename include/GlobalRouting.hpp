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

using namespace std ; 

int GlobalRoute(const vector<Bump>& bumps, const vector<double>& routing_area_coordinate, 
                const DesignRule& designRule, const string& output_path, vector<RoutingGraph>& allRDL, 
                double offset, vector<clock_t>& globalRouteTimes) ;

void ConstructRoutingGraph(const vector<Bump>& bumps, const vector<double> routing_area_coordinate, const DesignRule& designRule, const string& output_path, double Hor_SPACING_X, int layer, vector<RoutingGraph>& allRDL) ; 

vector<Bump> ProcessLeftAndRight(const vector<Bump>& bumps, const vector<Bump>& leftmostBumps, const vector<double>& routing_area_coordinate, double Hor_SPACING_X) ;
vector<Bump> FindLeftmostInEachRow(const vector<Bump>& bumps) ;
vector<EdgeNode> Triangulation(const vector<ViaNode>& via_nodes, const DesignRule& designRule) ;
