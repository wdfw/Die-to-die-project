#pragma once

// 定義 Via Node
struct ViaNode {
    string dieName;
    string type;
    int id;          
    double x, y;      

    ViaNode() : dieName(""), type(""), id(-1), x(0), y(0) {} 
    ViaNode(string dieName, string type, int id, double x, double y)
        : dieName(dieName), type(type), id(id), x(x), y(y) {}
};

// 定義 Edge Node，包含 Net-Sequence List (使用 Single Linked List)
struct EdgeNode {
    int id;                       
    double x, y;       
    vector<ViaNode> vias;      // 連接的兩個 ViaNode
    int capacity;                      // 邊的容量
    list<string> net_sequence;     // 單向鏈結串列 
    pair<double, double> start; // 存儲線段的起點
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



class RoutingGraph {
public:
    vector<ViaNode> via_nodes;              // 所有 Via Nodes
    vector<EdgeNode> edge_nodes;            // 所有 Edge Nodes
    vector<AccessViaEdge> access_via_edges; // 所有 Access-Via Edges
    vector<CrossTileEdge> cross_tile_edges; // 所有 Cross-Tile Edges
};
