#pragma once 

#include <vector>
#include <map>
#include <set>
#include <limits>
#include <random>


#include "DesignRule.hpp"
#include "Bump.hpp"
#include "Utils.hpp"
#include "RoutingGraph2.hpp"
#include "Router.hpp"
#include "Config.hpp"


extern vector<Bump> debugBumps ; 
extern vector<Net> debugNets ; 
extern vector<tuple<string, double, double>> debugLabels ; 
struct ChannelType ;
struct PhenotypeType ;
struct GenotypeType ; 
struct ChromosomeType ;

ostream& operator<<(ostream& os, const ChannelType& channel) ; 
ostream& operator<<(ostream& os, const PhenotypeType& phenotype) ; 
ostream& operator<<(ostream& os, const GenotypeType& genotype) ; 
ostream& operator<<(ostream& os, const ChromosomeType& chromosome) ; 

struct ChannelType{
    vector<TileToTileEdge> upChannel ;
    vector<TileToTileEdge> upCrossedChannel ;
    
    vector<TileToTileEdge> bottomChannel ;
    vector<TileToTileEdge> bottomCrossedChannel ;
};
    

struct PhenotypeType : public vector<TileToTileEdge> {
    int upChannelStartIndex ; 
    int upChannelEndIndex ;
    int bottomChannelStartIndex ;
    int bottomChannelEndIndex ;
    int channelNum ; 

    PhenotypeType(int upChannelStartIndex=-1, int upChannelEndIndex=-1, int bottomChannelStartIndex=-1, int bottomChannelEndIndex=-1, int channelNum=-1) : 
                upChannelStartIndex(upChannelStartIndex), upChannelEndIndex(upChannelEndIndex), bottomChannelStartIndex(bottomChannelStartIndex),
                bottomChannelEndIndex(bottomChannelEndIndex), channelNum(channelNum) {}

};

struct GenotypeType : public vector<int> {
    int startIndex ;
    int endIndex ;
    int codeLength ;
    GenotypeType(int startIndex=-1, int endIndex=-1, int codeLength=0) : startIndex(startIndex), endIndex(endIndex), codeLength(codeLength) {resize(codeLength);}
};

struct ChromosomeType{
    friend ostream& operator<<(ostream& os, const ChromosomeType& chromosome) ;
    void UpdatePhenotype(const vector<ChannelType>& channels) ;
    void Random(const vector<ChannelType>& channels) ;
    void Random() ; 
    
    PhenotypeType phenotype ;
    GenotypeType  genotype ;
    shared_ptr<vector<ChannelType>> channelsPtr ; 
};

struct ClusterChromosomes : vector<ChromosomeType>{
    double fitness = -1.0 ;  
};

class GARouter : public Router {


    double CacluteWireLength(ClusterChromosomes& chromosomes) ; 
    double CacluteConflictCount(ClusterChromosomes& chromosomes) ; 
    double CacluteCapacityValue(ClusterChromosomes& chromosomes) ; 
    double Fitness(ClusterChromosomes& chromosomes) ; 

    void Initial(const vector<tuple<int,int, int, int, int>>& genotypeInformations, const vector<tuple<int,int,int>>& phenotypeInformations, 
                 const vector<ChannelType>& channels, vector<ClusterChromosomes>& population) ; 
    
    void SelectParents(vector<ClusterChromosomes>& population, vector<ClusterChromosomes>& parents1, vector<ClusterChromosomes>& parents2) ; 
    void Crossover(vector<ClusterChromosomes>& population1, vector<ClusterChromosomes>& population2, vector<ClusterChromosomes>& offsprings) ; 
    void Mutation(vector<ClusterChromosomes>& population, vector<ClusterChromosomes>& offsprings) ; 
    void SelectSurviors(vector<ClusterChromosomes>& population,  vector<ClusterChromosomes>& offsprings, vector<ClusterChromosomes>& survivors) ;
    void Evaluate(vector<ClusterChromosomes>& population) ; 
    void EvaluateAll(vector<ClusterChromosomes>& population) ; 
    
    void ConstructChannel() ;  
    void ConstructRepresentation() ; 
public:
    default_random_engine generator ;
    bernoulli_distribution mutationDistrbution ;
    bernoulli_distribution crossoverDistrbution ;

    Configuration config ;
    using Router::Router ; 
    void GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph) override ;
} ;

// TileToTileEdge findRightUp(shared_ptr<ViaNode2> viaNode){
//     for(auto& node : viaNode->tileNodes){
//         if(node->y < viaNode->y && node->x > viaNode->x) return node ; 
//     }
//     return TileToTileEdge() ; 
// }

// TileToTileEdge findRightDown(shared_ptr<ViaNode2> viaNode){
//     for(auto& node : viaNode->tileNodes){
//         if(node->y > viaNode->y && node->x > viaNode->x) return node ; 
//     }
//     return TileToTileEdge() ; 
// }

// TileToTileEdge findLeftUp(shared_ptr<ViaNode2> viaNode){
//     for(auto& node : viaNode->tileNodes){
//         if(node->y < viaNode->y && node->x < viaNode->x) return node ; 
//     }
//     return TileToTileEdge() ; 
// }

// TileToTileEdge findLeftDown(shared_ptr<ViaNode2> viaNode){
//     for(auto& node : viaNode->tileNodes){
//         if(node->y > viaNode->y && node->x < viaNode->x) return node ; 
//     }
//     return TileToTileEdge() ; 
// }

// TileToTileEdge findLeft(shared_ptr<TileNode2> tileNode){
//     TileToTileEdge leftMostTileNode  ;

//     for(auto& node : tileNode->tileNodes){
//         if(node->x > tileNode->x){
//             if(leftMostTileNode==nullptr) leftMostTileNode = node ; 
//             else if(node->x > leftMostTileNode->x) leftMostTileNode = node ;  
//         }
//     }
//     return leftMostTileNode ; 
// }

// TileToTileEdge findUp(shared_ptr<TileNode2> tileNode){
//     TileToTileEdge upMostTileNode  ;

//     for(auto& node : tileNode->tileNodes){
//         if(node->y < tileNode->y){
//             if(upMostTileNode==nullptr) upMostTileNode = node ; 
//             else if(node->y < upMostTileNode->y) upMostTileNode = node ;  
//         }
//     }
//     return upMostTileNode ; 
// }

// class Router {


// void Router::GAChannelRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph){
//     vector<shared_ptr<ViaNode2>> routingViaNodes ; 
//     map<tuple<string, DieType, int>, shared_ptr<ViaNode2>> routingBumpsMapping ; 

//     map<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>> viaNodeMatching ; 
//     map<shared_ptr<ViaNode2>, vector<vector<TileToTileEdge>>> channels ; 

//     // map<shared_ptr<ViaNode2>, > viaNodeMatching ; 

//     for(auto& bump : routingBumps){
//         for(auto& viaNode : graph.viaNodes){
//             if(viaNode->name==bump.name && viaNode->type==bump.type && viaNode->id==bump.id && viaNode->type==SIGNAL){
//                 routingBumpsMapping[{bump.name, bump.type, bump.id}] = viaNode ; 
//             }
//         }
//     }

//     for(auto& [info1, viaNode1] : routingBumpsMapping){
//         for(auto& [info2, viaNode2] : routingBumpsMapping){
//             if(get<0>(info1)=="DIE1" && get<0>(info2)=="DIE2" && get<1>(info1)==get<1>(info2) && get<2>(info1)==get<2>(info2)){
//                 viaNodeMatching[viaNode1] = viaNode2 ; 
//             }
//         }
//     }

//     for(auto& [startViaNode, targetViaNode] : viaNodeMatching){
//         TileToTileEdge startTileNode = findRightUp(startViaNode) ;
//         TileToTileEdge targetTileNode = findLeftUp(targetViaNode) ;
//         TileToTileEdge currentTileNode = startTileNode ;
//         channels[startViaNode].resize(4) ; // 0, 0->1, 1, 1->0

//         channels[startViaNode][0].push_back(startTileNode) ; 
//         while(currentTileNode!=targetTileNode){
//             TileToTileEdge nextTileNode = findLeft(currentTileNode) ; 
//             // debugNets.push_back(Net("DUMMY", {{currentTileNode->x , currentTileNode->y, nextTileNode->x ,nextTileNode->y}})) ;
//             channels[startViaNode][0].push_back(nextTileNode) ; 
//             currentTileNode = nextTileNode ;
//         }
//         channels[startViaNode][1].resize(channels[startViaNode][0].size()) ;

//         startTileNode = findRightDown(startViaNode) ;
//         targetTileNode = findLeftDown(targetViaNode) ;
//         currentTileNode = startTileNode ;

//         channels[startViaNode][2].push_back(startTileNode) ; 
//         while(currentTileNode!=targetTileNode){
//             TileToTileEdge nextTileNode = findLeft(currentTileNode) ; 
//             // debugNets.push_back(Net("DUMMY", {{currentTileNode->x , currentTileNode->y, nextTileNode->x ,nextTileNode->y}})) ;
//             channels[startViaNode][2].push_back(nextTileNode) ; 
//             currentTileNode = nextTileNode ;
//         }
//         channels[startViaNode][3].resize(channels[startViaNode][2].size()) ;
        
//         map<TileToTileEdge, int> upChannelTileNodes ; 
//         for(int i=0; i<channels[startViaNode][0].size(); ++i){
//             upChannelTileNodes[channels[startViaNode][0][i]] = i ; 
//         }

//         for(int i=0; i<channels[startViaNode][2].size(); ++i){
//             for(auto& tileNode : channels[startViaNode][2][i]->tileNodes){
//                 if(upChannelTileNodes.find(tileNode)!=upChannelTileNodes.end()){
//                     channels[startViaNode][1][upChannelTileNodes[tileNode]] = channels[startViaNode][2][i] ;
//                     channels[startViaNode][3][i] = tileNode ;
//                     // debugNets.push_back(Net("DUMMY", {{tileNode->x , tileNode->y, channels[startViaNode][1][i]->x ,channels[startViaNode][1][i]->y}})) ;
//                 }
//             }
//         }

//         for(int i=0; i<channels[startViaNode][2].size()-1; ++i){
//             debugNets.push_back(Net("DUMMY", {{channels[startViaNode][2][i]->x , channels[startViaNode][2][i]->y, channels[startViaNode][2][i+1]->x , channels[startViaNode][2][i+1]->y}})) ;
//         }
//         for(int i=0; i<channels[startViaNode][2].size(); ++i){
//             if(channels[startViaNode][3][i])
//                 debugNets.push_back(Net("DUMMY", {{channels[startViaNode][2][i]->x , channels[startViaNode][2][i]->y, channels[startViaNode][3][i]->x , channels[startViaNode][3][i]->y}})) ;
//         }
//     }
    
//     //建立Vianode Representation的長度與意義
//     vector<vector<shared_ptr<ViaNode2>>> viaNodesMatix ; 
//     vector<int> represetSize ;
//     vector<vector<pair<int,int>>> represetSegement ;

//     vector<shared_ptr<ViaNode2>> sortedViaNodes ;
    
//     for(auto& [viaNode1, viaNode2] : viaNodeMatching) sortedViaNodes.push_back(viaNode1) ;
//     sort(sortedViaNodes.begin(), sortedViaNodes.end(), [](const shared_ptr<ViaNode2>& a, const shared_ptr<ViaNode2>& b) {return a->y < b->y;});

//     for (const auto& viaNode : sortedViaNodes) {
//         bool added = false;
//         for (auto& viaNodeRow : viaNodesMatix) {
//             if (!viaNodeRow.empty() && fabs(viaNodeRow[0]->y - viaNode->y) < epsilonY) {
//                 viaNodeRow.push_back(viaNode);
//                 added = true;
//                 break;
//             }
//         }
//         if (!added) viaNodesMatix.push_back({viaNode}) ;
//     }

//     represetSize = vector<int>(viaNodesMatix.size(), 0) ; 
//     represetSegement = vector<vector<pair<int,int>>>(viaNodesMatix.size()) ; 

//     for(int i=0, z; i<viaNodesMatix.size(); ++i){
//         set<TileToTileEdge> crossedEdge ; 
//         for(auto& viaNode : viaNodesMatix[i]){
//             for(auto& tileNode : channels[viaNode][1]) if(tileNode) crossedEdge.insert(tileNode) ; 
//         }
//         represetSize[i] = crossedEdge.size() ; 
//         represetSegement[i] = vector<pair<int,int>>(represetSize[i], {-1,-1}) ;

//         z = 0 ;
//         for(int j=0; j<viaNodesMatix[i].size(); ++j){
//             auto& viaNode = viaNodesMatix[i][j] ;
//             if(j!=0){
//                 auto& preViaNode = viaNodesMatix[i][j-1] ;
//                 for(int k=0; channels[preViaNode][0][k]!=findRightUp(viaNode); ++k){
//                     if(channels[preViaNode][1][k]) ++z ; 
//                 }
//             }
//             represetSegement[i][j].first = z ;
//             represetSegement[i][j].second = z-1 ;
//             for(auto& tileNode : channels[viaNode][1]){
//                 if(tileNode) ++represetSegement[i][j].second ;
//             }
//             // cout << represetSegement[i][j].first << " " << represetSegement[i][j].second << "\n" ;
//         }
//     }

//     //初始化Vianode的群集
//     struct GenotypeType{
//         shared_ptr<ViaNode2> viaNode ; 
//         int start ;
//         vector<int> numSequence ; 
//     };
//     struct PhenotypeType{
//         vector<TileToTileEdge> netSequence ; 
//     };
    

//     // vector<shared_ptr<ViaNode2>> viaNodesGroup ; for(auto& rowViaNodes : viaNodesMatix) viaNodesGroup.insert(viaNodesGroup.end(), rowViaNodes.begin(), rowViaNodes.end()) ;
//     vector<GenotypeType> genotypes ;
//     vector<GenotypeType> phenotypes ;

//     for(int i=0; i<viaNodesMatix.size(); ++i){
//         for(int j=0; j<viaNodesMatix[i].size(); ++j){
//             genotypes.push_back({viaNodesMatix[i][j], represetSegement[i][j].first, vector<int>(represetSize[i], 0)}) ;
//         }
//     }

// }