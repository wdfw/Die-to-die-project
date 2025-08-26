#pragma once
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DesignRule.hpp"
#include "Bump.hpp"
#include "Utils.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

using namespace std ; 

int RoutingGuideGenerate(vector<Bump> vias, vector<double> routing_area_coordinate, DesignRule designRule, string output_path, vector<RoutingGraph>& allRDL, double &offset, vector<int>& GlobalRouteTime) ;