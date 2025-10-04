#include "GAChannelRoute.hpp"
ostream& operator<<(ostream& os, const ChannelType& channel){
    for(auto& tileNode : channel.upChannel) os << tileNode->id << " " ;
    os << "\n" ;

    for(auto& tileNode : channel.upCrossedChannel){
        if(tileNode) os << tileNode->id << " " ;
        else os << "NULL " ; 
    }
    os << "\n" ;

    for(auto& tileNode : channel.bottomChannel) os << tileNode->id << " " ;
    os << "\n" ;

    for(auto& tileNode : channel.bottomCrossedChannel){
        if(tileNode) os << tileNode->id << " " ;
        else os << "NULL " ; 
    }
    os << "\n" ;
    return os ; 
}

ostream& operator<<(ostream& os, const PhenotypeType& phenotype){
    os << "(" << phenotype.size() << "|" << phenotype.upChannelStartIndex << "," << phenotype.upChannelEndIndex << "," <<  phenotype.bottomChannelStartIndex << "," << phenotype.bottomChannelEndIndex << "," << phenotype.channelNum  << ")" ;
    for(auto& v : phenotype){
        if(v) os << v->id << " "  ;
        else os << "NULL " ;
    }
    for(int i=0; i<phenotype.size()-1; ++i){
         debugNets.push_back(Net("ground", {{phenotype[i]->x ,phenotype[i]->y, phenotype[i+1]->x, phenotype[i+1]->y}})) ;
    }
    return os ; 
}   

ostream& operator<<(ostream& os, const GenotypeType& genotype){
    os << "(" << genotype.startIndex << "," << genotype.endIndex << "," <<  genotype.codeLength  << ")" ;
    for(auto& v : genotype) os << v ;
    return os ; 
}

ostream& operator<<(ostream& os, const ChromosomeType& chromosome){
    os << chromosome.phenotype << "\n" ;
    os << chromosome.genotype << "\n" ;
    return os ; 
}



void ChromosomeType::UpdatePhenotype(const vector<ChannelType>& channels){
    int s=-1, u=phenotype.upChannelStartIndex, d=phenotype.bottomChannelStartIndex ;
    const ChannelType& channel = channels[phenotype.channelNum] ;

    for(int i=genotype.startIndex; u<phenotype.upChannelEndIndex && d<phenotype.bottomChannelEndIndex && i<genotype.endIndex; ){
        if(s==-1){
            if(genotype[i]==0){
                phenotype.push_back(channel.upChannel[u++]) ;
                s = 0 ; 
            }else if(genotype[i]==1){
                phenotype.push_back(channel.bottomChannel[d++]) ;
                s = 1 ;
            }
            ++i ;
        }else{
            if(s==0){
                phenotype.push_back(channel.upChannel[u]) ;
                if(channel.upCrossedChannel[u]){
                    if(genotype[i]==1){
                        while(channel.bottomChannel[d]!=channel.upCrossedChannel[u]){
                            ++d ; 
                        }
                        phenotype.push_back(channel.bottomChannel[d]) ;
                        ++d ; 
                        s = 1 ; 
                    }
                    ++i ;
                }
                ++u ; 
            }else if(s==1){
                phenotype.push_back(channel.bottomChannel[d]) ;
                if(channel.bottomCrossedChannel[d]){
                    if(genotype[i]==0){
                        while(channel.upChannel[u]!=channel.bottomCrossedChannel[d]){
                            ++u ; 
                        }
                        phenotype.push_back(channel.upChannel[u]) ;
                        ++u ; 
                        s = 0 ; 
                    }
                    ++i ;
                }
                ++d ; 
            }
        }

    }
}

void ChromosomeType::Random(const vector<ChannelType>& channels){
    Random() ; 
    UpdatePhenotype(channels) ; 
}

void ChromosomeType::Random(){
    for(int i=0; i<genotype.size(); ++i){
        if(genotype.startIndex<=i && i<genotype.endIndex) genotype[i] = rand()%2 ; 
        else genotype[i] = 2 ; 
    }
}

TileToTileEdge findRightUp(shared_ptr<ViaNode2> viaNode){
    for(auto& node : viaNode->tileNodes){
        if(node->y < viaNode->y && node->x > viaNode->x) return node ; 
    }
    return TileToTileEdge() ; 
}

TileToTileEdge findRightDown(shared_ptr<ViaNode2> viaNode){
    for(auto& node : viaNode->tileNodes){
        if(node->y > viaNode->y && node->x > viaNode->x) return node ; 
    }
    return TileToTileEdge() ; 
}

TileToTileEdge findLeftUp(shared_ptr<ViaNode2> viaNode){
    TileToTileEdge leftMostTileNode ; 
    for(auto& node : viaNode->tileNodes){
        if(node->y < viaNode->y && node->x < viaNode->x){
            if(!leftMostTileNode) leftMostTileNode = node ; 
            else if(leftMostTileNode->x < viaNode->x) leftMostTileNode = node ; 
        }
    }
    return leftMostTileNode ; 
}

TileToTileEdge findLeftDown(shared_ptr<ViaNode2> viaNode){
    TileToTileEdge leftMostTileNode ; 
    for(auto& node : viaNode->tileNodes){
        if(node->y > viaNode->y && node->x < viaNode->x){
            if(!leftMostTileNode) leftMostTileNode = node ; 
            else if(leftMostTileNode->x < viaNode->x) leftMostTileNode = node ; 
        }
    }
    return leftMostTileNode ; 
}

TileToTileEdge findLeft(shared_ptr<TileNode2> tileNode){
    TileToTileEdge leftMostTileNode  ;

    for(auto& node : tileNode->tileNodes){
        if(node->x > tileNode->x){
            if(leftMostTileNode==nullptr) leftMostTileNode = node ; 
            else if(node->x > leftMostTileNode->x) leftMostTileNode = node ;  
        }
    }
    return leftMostTileNode ; 
}

TileToTileEdge findUp(shared_ptr<TileNode2> tileNode){
    TileToTileEdge upMostTileNode  ;

    for(auto& node : tileNode->tileNodes){
        if(node->y < tileNode->y){
            if(upMostTileNode==nullptr) upMostTileNode = node ; 
            else if(node->y < upMostTileNode->y) upMostTileNode = node ;  
        }
    }
    return upMostTileNode ; 
}


void GARouter::GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph){
    map<tuple<string, DieType, int>, shared_ptr<ViaNode2>> routingBumpsMapping ; 
    map<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>> viaNodeMatching ; 

    for(auto& bump : routingBumps){
        for(auto& viaNode : graph.viaNodes){
            if(viaNode->name==bump.name && viaNode->type==bump.type && viaNode->id==bump.id && viaNode->type==SIGNAL){
                routingBumpsMapping[{bump.name, bump.type, bump.id}] = viaNode ; 
            }
        }
    }

    for(auto& [info1, viaNode1] : routingBumpsMapping){
        for(auto& [info2, viaNode2] : routingBumpsMapping){
            if(get<0>(info1)=="DIE1" && get<0>(info2)=="DIE2" && get<1>(info1)==get<1>(info2) && get<2>(info1)==get<2>(info2)){
                viaNodeMatching[viaNode1] = viaNode2 ; 
            }
        }
    }

    vector<vector<shared_ptr<ViaNode2>>> die1ViaNodesMatix, die2ViaNodesMatix, combinedViaNodesMatix ; 
    vector<shared_ptr<ViaNode2>> die1RoutingViaNodes, die2RoutingViaNodes ; 
    vector<int> representSize ; 
    vector<tuple<int,int,int>> genotypeInformations ; 
    vector<tuple<int,int, int, int, int>> phenotypeInformations ; 

    vector<ChannelType> channels ; 

    for(auto& [viaNode1, viaNode2] : viaNodeMatching){
        die1RoutingViaNodes.push_back(viaNode1) ;
        die2RoutingViaNodes.push_back(viaNode2) ;
    }

    Matrixization(die1RoutingViaNodes, die1ViaNodesMatix, epsilonY) ;
    Matrixization(die2RoutingViaNodes, die2ViaNodesMatix, epsilonY) ;

    if(die1RoutingViaNodes.size() != die2RoutingViaNodes.size()) cerr << "Unmatched routing vias." ;
    for(int i=0; i<die2ViaNodesMatix.size(); ++i){
        combinedViaNodesMatix.push_back(die1ViaNodesMatix[i]) ; 
        combinedViaNodesMatix.back().insert(combinedViaNodesMatix.back().end(), die2ViaNodesMatix[i].begin(), die2ViaNodesMatix[i].end()) ;
    }
    for(int i=0, crossedCount; i<combinedViaNodesMatix.size(); ++i){
        shared_ptr<ViaNode2> leftmostViaNode = combinedViaNodesMatix[i].front() ; 
        shared_ptr<ViaNode2> rightmostViaNode = combinedViaNodesMatix[i].back() ; 
        map<TileToTileEdge, int> upChannelTileNodeIndexes ; 
        
        TileToTileEdge startTileNode = findRightUp(leftmostViaNode) ;
        TileToTileEdge targetTileNode = findLeftUp(rightmostViaNode) ;
        TileToTileEdge currentTileNode = startTileNode ;
        channels.push_back({}) ; 

        channels.back().upChannel.push_back(currentTileNode) ; 
        while(currentTileNode!=targetTileNode){
            TileToTileEdge nextTileNode = findLeft(currentTileNode) ; 
            channels.back().upChannel.push_back(nextTileNode) ; 
            currentTileNode = nextTileNode ;
        }        
        
        startTileNode = findRightDown(leftmostViaNode) ;
        targetTileNode = findLeftDown(rightmostViaNode) ;
        currentTileNode = startTileNode ;

        channels.back().bottomChannel.push_back(currentTileNode) ; 
        while(currentTileNode!=targetTileNode){
            TileToTileEdge nextTileNode = findLeft(currentTileNode) ; 
            channels.back().bottomChannel.push_back(nextTileNode) ; 
            currentTileNode = nextTileNode ;
        }

        for(int i=0; i<channels.back().upChannel.size(); ++i) upChannelTileNodeIndexes[channels.back().upChannel[i]] = i ; 
        
        channels.back().upCrossedChannel = vector<TileToTileEdge>(channels.back().upChannel.size()) ;
        channels.back().bottomCrossedChannel = vector<TileToTileEdge>(channels.back().bottomChannel.size()) ;
        crossedCount = 0 ;
        for(int i=0; i<channels.back().bottomChannel.size(); ++i){
            for(auto& tileNode : channels.back().bottomChannel[i]->tileNodes){
                if(upChannelTileNodeIndexes.find(tileNode)!=upChannelTileNodeIndexes.end()){
                    channels.back().bottomCrossedChannel[i] = tileNode ; 
                    channels.back().upCrossedChannel[upChannelTileNodeIndexes[tileNode]] = channels.back().bottomChannel[i] ;
                    ++crossedCount ;
                }
            }
        }

        representSize.push_back(crossedCount) ; 

        // for(int i=0; i<channels.back().upChannel.size()-1; ++i){
            // debugNets.push_back(Net("DUMMY", {{channels.back().upChannel[i]->x , channels.back().upChannel[i]->y, channels.back().upChannel[i+1]->x , channels.back().upChannel[i+1]->y}})) ;
        // }
        // for(int i=0; i<channels.back().bottomChannel.size()-1; ++i){
        //     debugNets.push_back(Net("DUMMY", {{channels.back().bottomChannel[i]->x , channels.back().bottomChannel[i]->y, channels.back().bottomChannel[i+1]->x , channels.back().bottomChannel[i+1]->y}})) ;
        // }
    }
    

    for(int i=0; i<die1ViaNodesMatix.size(); ++i){
        for(int j=0; j<die1ViaNodesMatix[i].size(); ++j){
            int startIndex = -1, targetIndex = -1, currentIndex = 0 ;
            int upChannelStartIndex = -1, bottomChannelStartIndex = -1, upChannelEndIndex = -1, bottomChannelEndIndex = -1, channelNum = i ;
            shared_ptr<ViaNode2> stratViaNode = die1ViaNodesMatix[i][j] ; 
            shared_ptr<ViaNode2> targetViaNode = viaNodeMatching[stratViaNode] ; 
            TileToTileEdge startTileNode = findRightUp(stratViaNode) ;
            TileToTileEdge targetTileNode = findLeftUp(targetViaNode) ;
            for(int k=0; k<channels[i].upChannel.size(); ++k){
                if(channels[i].upChannel[k]==startTileNode){
                    startIndex = currentIndex ; 
                }else if(channels[i].upChannel[k]==targetTileNode){
                    targetIndex = currentIndex+1 ; 
                    break ;
                }
                if(channels[i].upCrossedChannel[k]) ++currentIndex ;
            }

            startTileNode = findRightUp(stratViaNode) ;
            targetTileNode = findLeftUp(targetViaNode) ;
            for(int k=0; k<channels[i].upChannel.size(); ++k){
                if(channels[i].upChannel[k]==startTileNode){
                    upChannelStartIndex = k ; 
                }else if(channels[i].upChannel[k]==targetTileNode){
                    upChannelEndIndex = k+1 ; 
                    break ;
                }
            }

            startTileNode = findRightDown(stratViaNode) ;
            targetTileNode = findLeftDown(targetViaNode) ;
            for(int k=0; k<channels[i].bottomChannel.size(); ++k){
                if(channels[i].bottomChannel[k]==startTileNode){
                    bottomChannelStartIndex = k ; 
                }else if(channels[i].bottomChannel[k]==targetTileNode){
                    bottomChannelEndIndex = k+1 ; 
                    break ;
                }
            }
            genotypeInformations.push_back({startIndex, targetIndex, representSize[i]}) ;
            phenotypeInformations.push_back({upChannelStartIndex, upChannelEndIndex, bottomChannelStartIndex, bottomChannelEndIndex, channelNum}) ;
        }
    }
    //------------------------------------------------------------------------------------------------
    //-------------------------------------GA Start---------------------------------------------------
    //------------------------------------------------------------------------------------------------
    vector<ClusterChromosomes> population ;
    map<int, shared_ptr<ViaNode2>> resultIndexMapping ; 
    
    //Initialize
    Initial(phenotypeInformations, genotypeInformations, channels, population) ;
}
void GARouter::Initial(const vector<tuple<int,int, int, int, int>>& phenotypeInformations, const vector<tuple<int,int,int>>& genotypeInformations, const vector<ChannelType>& channels, vector<ClusterChromosomes>& population){
    populations = vector<ClusterChromosomes>(config.population_size) ; 

    for(int i=0; i<populations.size(); ++i){
        ClusterChromosomes& cluster = populations[i] ;
       
        for(int j=0; j<phenotypeInformations.size(); ++j){
            ChromosomeType chromosome = {{
                get<0>(phenotypeInformations[j]),
                get<1>(phenotypeInformations[j]),
                get<2>(phenotypeInformations[j]),
                get<3>(phenotypeInformations[j]),
                get<4>(phenotypeInformations[j])
            }, {
                get<0>(genotypeInformations[j]),
                get<1>(genotypeInformations[j]),
                get<2>(genotypeInformations[j])
            }} ;

            chromosome.channelsPtr = make_shared<vector<ChannelType>>(channels) ; 
            chromosome.Random(*chromosome.channelsPtr) ; 
            cluster.push_back(chromosome) ; 
        }
    }
}

double GARouter::Fitness(ClusterChromosomes& chromosomes){
    
}

void GARouter::Evalute(vector<ClusterChromosomes>& population){
    for(auto& cluster : population){
        if(cluster.fitness<0.0) cluster.fitness = Fitness(cluster) ; 
    }
}
