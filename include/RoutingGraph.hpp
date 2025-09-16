#pragma once
#include "Bump.hpp"


// 定義 Via Node
struct ViaNode {
    string dieName;
    DieType type;
    int id;          
    double x, y;      
    ViaNode(const string& dieName="", DieType type=DUMMY, int id=-1, double x=0.0, double y=0.0) : dieName(dieName), type(type), id(id), x(x), y(y) {}
};

// 定義 Edge Node，包含 Net-Sequence List (使用 Single Linked List)
struct EdgeNode {
    int id;                       
    double x, y;       
    vector<ViaNode> vias;      // 連接的兩個 ViaNode
    int capacity;                      // 邊的容量
    list<string> net_sequence;     // 單向鏈結串列 
    pair<double, double> start; // 存儲線段的起點z
    pair<double, double> end;  // 存儲線段的終點

    EdgeNode(int id, double x, double y, vector<ViaNode> vias, int capacity)
        : id(id), x(x), y(y), vias(vias), capacity(capacity) {}
};

// 定義 Access-Via Edge
struct AccessViaEdge {
    int id;          
    ViaNode via;     // 連接的 ViaNode
    EdgeNode edge;     // 連接的 EdgeNode

    AccessViaEdge(int id, ViaNode via, EdgeNode edge)
        : id(id), via(via), edge(edge) {}
};

// 定義 Cross-Tile Edge，包含角度資訊
struct CrossTileEdge {
    int id;                   
    vector<EdgeNode> edges;  // 連接的兩個 EdgeNode
    double angle;             // 角度

    CrossTileEdge(int id, vector<EdgeNode> edges, double angle)
        : id(id), edges(edges), angle(angle) {}
};


struct Edge{
    int id;                       
    vector<ViaNode> vias;      // 連接的兩個 ViaNode
    int capacity;                      // 邊的容量
    pair<double, double> start; // 存儲線段的起點z
    pair<double, double> end;  // 存儲線段的終點

    Edge(int id, vector<ViaNode> vias, int capacity)
        : id(id), vias(vias), capacity(capacity) {
            start = {vias[0].x, vias[0].y} ;
            end = {vias[1].x, vias[1].y} ;
        }
};

struct TruthEdgeNode {
    int id;                       
    double x, y;       
    vector<ViaNode> vias;      // 連接的三個 ViaNode
    vector<Edge> edges;      // 連接的三個 ViaNode

    TruthEdgeNode(int id, double x, double y, const vector<ViaNode>& vias, const vector<Edge>& edges): id(id), x(x), y(y), vias(vias), edges(edges) {}
};


class RoutingGraph {
public:
    vector<ViaNode> via_nodes;              // 所有 Via Nodes
    
    vector<TruthEdgeNode> _edge_nodes;            // 所有 Edge Nodes
    vector<EdgeNode> edge_nodes;            // 所有 Edge Nodes
    
    vector<AccessViaEdge> access_via_edges; // 所有 Access-Via Edges // via to tile
    vector<CrossTileEdge> cross_tile_edges; // 所有 Cross-Tile Edges // tile to tile
};

