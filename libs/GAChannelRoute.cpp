#include "GAChannelRoute.hpp"
bool showExpand = false ; 
bool ssh=false ; 
const int ReportInterval = 100 ;

ostream& operator<<(ostream& os, const ChannelType& channel){
    int i = 0 ;
    for(auto& tileNode : channel.upChannel) os << tileNode->id << "," << i++ << " " ;
    os << "\n" ;
    i=0 ;
    for(auto& tileNode : channel.upCrossedChannel){
        if(tileNode) os << tileNode->id << "," << i << " " ;
        else os << "NULL " ; 
        ++i ; 

    }
    os << "\n" ;
    i=0;
    for(auto& tileNode : channel.bottomChannel) os << tileNode->id << "," << i++ << " " ;
    os << "\n" ;
    i=0 ; 
    for(auto& tileNode : channel.bottomCrossedChannel){
        if(tileNode) os << tileNode->id << "," << i << " " ;
        else os << "NULL " ; 
        ++i ; 
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
    // for(int i=0; i<phenotype.size()-1; ++i){
    //      debugNets.push_back(Net("ground", {{phenotype[i]->x ,phenotype[i]->y, phenotype[i+1]->x, phenotype[i+1]->y}})) ;
    // }
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


void ChromosomeType::LegalizeGenotype(){
    genotype[genotype.endIndex-1] = genotype[genotype.endIndex-2] ;
}



void ChromosomeType::UpdatePhenotype(const vector<ChannelType>& channels){
    int s=-1, u=phenotype.upChannelStartIndex, d=phenotype.bottomChannelStartIndex ;
    // cout << "--------------------------------------\n" ;
    // for(int i=0; i<channels.size(); ++i) cout << channels[i] << "\n" ;
    // throw runtime_error("") ;
    const ChannelType& channel = channels[phenotype.channelNum] ;
    phenotype.clear() ;

                // for(auto& tileNode : channel.upChannel) os << tileNode->id << "," << i++ << " " ;
                // for(auto& tileNode : channel.upCrossedChannel){
                // for(auto& tileNode : channel.bottomChannel) os << tileNode->id << "," << i++ << " " ;
                // for(auto& tileNode : channel.bottomCrossedChannel){
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
                        phenotype.push_back(channel.upCrossedChannel[u]) ;
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
                        phenotype.push_back(channel.bottomCrossedChannel[d]) ;
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
    LegalizeGenotype() ; 
}

void ChromosomeType::Random(){
    for(int i=0; i<genotype.size(); ++i){
        if(genotype.startIndex<=i && i<genotype.endIndex) genotype[i] = rand()%2 ; 
        else genotype[i] = 2 ; 
    }
}

void ClusterChromosomes::ToGlobalNet(const vector<int>& orders, vector<GlobalNet>& nets){
    nets.clear() ;
    for(int i=0,j; i<orders.size(); ++i){
        j = orders[i] ;
        if(at(j).startViaNode){
            shared_ptr<ViaNode2> viaNode = at(j).startViaNode ;
            nets.push_back(GlobalNet{viaNode->name, viaNode->type, viaNode->id, viaNode, at(j).endViaNode}) ;
        }else{
            nets.push_back(GlobalNet{"", VSS, -1}) ;
        }
        
        // cout << *at(j).startViaNode << ":\n" ;
        for(int k=0; k<at(j).phenotype.size(); ++k){
            // cout << at(j).phenotype[k] << " " << at(j).phenotype[k]->id << " " ;
            // if(at(j).phenotype[k].crossedViaNode1) cout << *at(j).phenotype[k].crossedViaNode1 << " " << *at(j).phenotype[k].crossedViaNode2 << "|" ;
            nets.back().push_back(at(j).phenotype[k]) ;
        }
        // cout << "\n" ;
    }
}


TileToTileEdge findRightUp(shared_ptr<ViaNode2> viaNode){
    double theta, closestTheta = 10.0 ;
    TileToTileEdge closestViaNode = TileToTileEdge() ;
    for(auto& node : viaNode->tileNodes){
        if(node->y < viaNode->y && node->x > viaNode->x){
            theta = atan2(viaNode->y - node->y, node->x - viaNode->x) ;
            if( fabs(theta - M_PI/6) < fabs(closestTheta) ){ // -30度
                closestTheta = theta - M_PI/6 ;  
                closestViaNode = node ;
            }
        }
    }
    return closestViaNode ; 
}

TileToTileEdge findRightDown(shared_ptr<ViaNode2> viaNode){
    double theta, closestTheta = 10.0 ;
    TileToTileEdge closestViaNode = TileToTileEdge() ;
    for(auto& node : viaNode->tileNodes){
        if(node->y > viaNode->y && node->x > viaNode->x){
            theta = atan2(node->y - viaNode->y, node->x - viaNode->x) ;
            if( fabs(theta - M_PI/6) < fabs(closestTheta) ){ // -30度
                closestTheta = theta - M_PI/6 ;  
                closestViaNode = node ;
            }
        }
    }
    return closestViaNode ; 
}

TileToTileEdge findLeftUp(shared_ptr<ViaNode2> viaNode){
    double theta, closestTheta = 10.0 ;
    TileToTileEdge closestViaNode = TileToTileEdge() ;
    for(auto& node : viaNode->tileNodes){
        if(node->y < viaNode->y && node->x < viaNode->x){
            theta = atan2(viaNode->y - node->y, viaNode->x - node->x) ;
            if( fabs(theta - M_PI/6) < fabs(closestTheta) ){ // 120度
                closestTheta = theta - M_PI/6 ;  
                closestViaNode = node ;
            }
        }
    }
    return closestViaNode ; 
}

TileToTileEdge findLeftDown(shared_ptr<ViaNode2> viaNode){
    double theta, closestTheta = 10.0 ;
    TileToTileEdge closestViaNode = TileToTileEdge() ;
    for(auto& node : viaNode->tileNodes){
        if(node->y > viaNode->y && node->x < viaNode->x){
            theta = atan2(node->y - viaNode->y, viaNode->x - node->x) ;
            if( fabs(theta - M_PI/6) < fabs(closestTheta) ){ // 270度
                closestTheta = theta - M_PI/6 ;  
                closestViaNode = node ;
            }
        }
    }
    return closestViaNode ; 
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

void GARouter::AddShielding(ClusterChromosomes& chromosomes, vector<int>& orders){
    int crrentSize = chromosomes.size() ; 
    vector<int> newOrders ; 
    for(int i=0, j=0, k=0; i<orders.size(); ++i){
        j = orders[i] ; 
        ChromosomeType upperChromosome = chromosomes[j] ; 
        ChromosomeType bottomChromosome = chromosomes[j] ; 
        
        upperChromosome.startViaNode = nullptr ; 
        upperChromosome.endViaNode = nullptr ; 
        bottomChromosome.startViaNode = nullptr ; 
        bottomChromosome.endViaNode = nullptr ; 


        chromosomes.push_back(upperChromosome) ; newOrders.push_back(chromosomes.size()-1) ; 
        newOrders.push_back(j) ; 
        chromosomes.push_back(bottomChromosome) ; newOrders.push_back(chromosomes.size()-1) ; 
    }

    orders = newOrders ; 
}

double GARouter::GlobalRoute(const vector<Bump>& routingBumps, RoutingGraph2& graph, vector<GlobalNet>& nets){
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
    vector<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>> connections ; 

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
        map<TileToTileEdge, int> bottomChannelTileNodeIndexes ; 
        
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
        
        map<TileToTileEdge, TileToTileEdge> crossedTileToTileEdge ;

        for(auto edge : channels.back().upChannel){
            crossedTileToTileEdge[edge] = edge ; 
            for(auto& tileNode : graph.tileNodes){
                for(auto edge2 : tileNode->tileNodes){
                    if( edge==edge2 && edge.crossedViaNode1 && edge.crossedViaNode2 && edge2.crossedViaNode1 && edge2.crossedViaNode2 &&
                        !((edge.crossedViaNode1==edge2.crossedViaNode1 && edge.crossedViaNode2==edge2.crossedViaNode2) ||
                         (edge.crossedViaNode1==edge2.crossedViaNode2 && edge.crossedViaNode2==edge2.crossedViaNode2))
                    ){
                        if(fabs(edge2.crossedViaNode2->y - edge2.crossedViaNode1->y) < 1.0) crossedTileToTileEdge[edge] = edge2;
                    }
                }
            }
        }

        for(auto edge : channels.back().bottomChannel){
            crossedTileToTileEdge[edge] = edge ; 
            for(auto& tileNode : graph.tileNodes){
                for(auto edge2 : tileNode->tileNodes){
                    if( edge==edge2 && edge.crossedViaNode1 && edge.crossedViaNode2 && edge2.crossedViaNode1 && edge2.crossedViaNode2 &&
                        !((edge.crossedViaNode1==edge2.crossedViaNode1 && edge.crossedViaNode2==edge2.crossedViaNode2) ||
                         (edge.crossedViaNode1==edge2.crossedViaNode2 && edge.crossedViaNode2==edge2.crossedViaNode2))
                    ){
                        if(fabs(edge2.crossedViaNode2->y - edge2.crossedViaNode1->y) < 1.0) crossedTileToTileEdge[edge] = edge2;

                    }
                }
            }
        }
        
        // for(auto& p : crossedTileToTileEdge){
        //     cout << p.first << " " << p.first->id << " "  ;
        //     if(p.first.crossedViaNode1)cout << *p.first.crossedViaNode1 << " " << *p.first.crossedViaNode2 << "\n" ;
        //     else cout << "\n" ;
        //     cout << p.second << " " << p.second->id << " " ;
        //     if(p.second.crossedViaNode2 ) cout << *p.second.crossedViaNode1 << " " << *p.second.crossedViaNode2 << " " << *p.second.capacity << " pp\n" ;
        //     else cout << "\n" ;
            
        // }
        for(int i=0; i<channels.back().upChannel.size(); ++i) upChannelTileNodeIndexes[channels.back().upChannel[i]] = i ; 
        for(int i=0; i<channels.back().bottomChannel.size(); ++i) bottomChannelTileNodeIndexes[channels.back().bottomChannel[i]] = i ; 

        channels.back().upCrossedChannel = vector<TileToTileEdge>(channels.back().upChannel.size()) ;
        channels.back().bottomCrossedChannel = vector<TileToTileEdge>(channels.back().bottomChannel.size()) ;
        crossedCount = 0 ;
        
        for(int i=0; i<channels.back().bottomChannel.size(); ++i){
            for(auto& tileNode : channels.back().bottomChannel[i]->tileNodes){
                if(upChannelTileNodeIndexes.find(tileNode)!=upChannelTileNodeIndexes.end()){
                    channels.back().bottomCrossedChannel[i] = crossedTileToTileEdge[tileNode] ; 
                    channels.back().upCrossedChannel[upChannelTileNodeIndexes[tileNode]] = crossedTileToTileEdge[channels.back().bottomChannel[i]] ;
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

    // for(auto& cc:channels){
    //     for(auto& w : cc.upCrossedChannel) if(w) cout << w->id << " " << *w.crossedViaNode1 << " " << *w.crossedViaNode2 << "|\n" ;
    // }

    for(int i=0; i<die1ViaNodesMatix.size(); ++i){

        for(int j=0; j<die1ViaNodesMatix[i].size(); ++j){
            int startIndex = -1, targetIndex = -1, currentIndex = 0 ;
            int upChannelStartIndex = -1, bottomChannelStartIndex = -1, upChannelEndIndex = -1, bottomChannelEndIndex = -1, channelNum = i ;
            shared_ptr<ViaNode2> stratViaNode = die1ViaNodesMatix[i][j] ; 
            shared_ptr<ViaNode2> targetViaNode = viaNodeMatching[stratViaNode] ; 
            TileToTileEdge startTileNode = findRightUp(stratViaNode) ;
            TileToTileEdge targetTileNode = findLeftUp(targetViaNode) ;

            // cout << *stratViaNode << " " << *startTileNode << "|\n" ;

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

            connections.push_back({stratViaNode, targetViaNode}) ;
            genotypeInformations.push_back({startIndex, targetIndex, representSize[i]}) ;
            phenotypeInformations.push_back({upChannelStartIndex, upChannelEndIndex, bottomChannelStartIndex, bottomChannelEndIndex, channelNum}) ;
        }
    }
    //------------------------------------------------------------------------------------------------
    //-------------------------------------GA Start---------------------------------------------------
    //------------------------------------------------------------------------------------------------
    vector<shared_ptr<ClusterChromosomes>> population, parents1, parents2, crossoverOffsprings, mutationOffsprings, surviors ;
    map<int, shared_ptr<ViaNode2>> resultIndexMapping ; 
    
    //Initialize
    // Initial(phenotypeInformations, genotypeInformations, connections, channels, population) ;
    Initial(phenotypeInformations, genotypeInformations, connections, channels, population) ;

    
    // cout << "IT 0: " ; EvaluateAll(population) ;

// Execute Time Analysis:
// SelectParents: 0 ms
// Crossover: 1207 ms
// Mutation: 5577 ms
// UpdatePhenotype: 1339 ms
// Evaluate: 6240 ms
// SelectSurviors: 4 ms
// EvaluateAll: 5651 ms

// Execute Time Analysis:
// SelectParents: 3376 ms
// Crossover: 2203 ms
// Mutation: 1076 ms
// UpdatePhenotype: 1045 ms
// Evaluate: 6082 ms
// SelectSurviors: 3684 ms
// EvaluateAll: 5472 ms
    Timer timer ; timer.SetClock() ; 
     for(int i=0; i<config.num_generations ; i++){

        globalTimers[0].SetClock() ; 
        SelectParents(population, parents1, parents2) ;
        globalExecTimes[0] += globalTimers[0].GetDurationMilliseconds() ; 

        globalTimers[0].SetClock() ; 
        Crossover(parents1, parents2, crossoverOffsprings) ;
        globalExecTimes[1] += globalTimers[0].GetDurationMilliseconds() ; 

        globalTimers[0].SetClock() ; 
        Mutation(crossoverOffsprings, mutationOffsprings) ;
        globalExecTimes[2] += globalTimers[0].GetDurationMilliseconds() ; 
       
        globalTimers[0].SetClock() ; 
        for(auto& cc : mutationOffsprings){
            if(cc->fitness<0){
                for(auto& c : *cc) c.LegalizeGenotype() ;
                for(auto& c : *cc) c.UpdatePhenotype(*c.channelsPtr) ;
            }
        }
        globalExecTimes[3] += globalTimers[0].GetDurationMilliseconds() ; 

        globalTimers[0].SetClock() ; 
        Evaluate(mutationOffsprings) ; 
        globalExecTimes[4] += globalTimers[0].GetDurationMilliseconds() ; 

        globalTimers[0].SetClock() ; 
        SelectSurviors(population, mutationOffsprings, surviors) ; population = surviors ;
        globalExecTimes[5] += globalTimers[0].GetDurationMilliseconds() ; 

        globalTimers[0].SetClock() ; 
        EvaluateAll(population) ; 
        globalExecTimes[6] += globalTimers[0].GetDurationMilliseconds() ; 

        if((i+1)% ReportInterval==0){
            cout << "IT " + to_string(i+1) + ": Objective " ; 
            double sum = 0.0 ; 
            for(auto& c : population) sum += c->fitness ;
            cout << sum/population.size() << " | Elapse " << timer.GetDurationMilliseconds() << " ms" << " | Average " << timer.GetDurationMilliseconds()/50 << " ms\n" ;
            timer.SetClock() ; 
        }
        // ClusterChromosomes bestClusterChromosomes = *min_element(population.begin(), population.end(), [](ClusterChromosomes& c1, ClusterChromosomes& c2){return c1.fitness<c2.fitness;}) ;
        // if(bestClusterChromosomes.fitness==0) break;
    }


    cout << "Execute Time Analysis:" << "\n" ;
    cout << "SelectParents: " << globalExecTimes[0] << " ms" "\n" ;
    cout << "Crossover: " << globalExecTimes[1] << " ms" "\n" ;
    cout << "Mutation: " << globalExecTimes[2] << " ms" "\n" ;
    cout << "UpdatePhenotype: " << globalExecTimes[3] << " ms" "\n" ;
    cout << "Evaluate: " << globalExecTimes[4] << " ms" "\n" ;
    // cout << "   Conflict: " << globalExecTimes[7] << " ms" "\n" ;
    // cout << "   Capacity: " << globalExecTimes[8] << " ms" "\n" ;
    // cout << "   Wire length: " << globalExecTimes[9] << " ms" "\n" ;
    cout << "SelectSurviors: " << globalExecTimes[5] << " ms" "\n" ;
    cout << "EvaluateAll: " << globalExecTimes[6] << " ms" "\n\n" ;
    ClusterChromosomes& bestClusterChromosomes = **min_element(population.begin(), population.end(), 
                [](shared_ptr<ClusterChromosomes>& c1, shared_ptr<ClusterChromosomes>& c2){return c1->fitness<c2->fitness;}) ;

    // showExpand = true ;
   
    int bestConflictCount = CacluteConflictCount(bestClusterChromosomes) ;

    cout << "Best result: " << bestClusterChromosomes.fitness << "\n" ;
    cout << "Conflic count: " << bestConflictCount << "\n" ;
    cout << "Exceed capacity: " << CacluteCapacityValue(bestClusterChromosomes) << "\n" ;

    double feedback = -CacluteConflictCount(bestClusterChromosomes) ; 
    if(feedback==0.0 && CacluteCapacityValue(bestClusterChromosomes)) feedback = -1 ; 

    if(bestConflictCount==0){
        vector<int> orders ; 
        DetermineOrder(bestClusterChromosomes, orders) ;
        AddShielding(bestClusterChromosomes, orders) ; 
        bestClusterChromosomes.ToGlobalNet(orders, nets) ; 
    }

    // for(auto& c : bestClusterChromosomes) cout << c << "\n" ;
    return feedback ; 
}

void GARouter::SelectSurviors(vector<shared_ptr<ClusterChromosomes>>& population, 
                        vector<shared_ptr<ClusterChromosomes>>& offsprings,
                        vector<shared_ptr<ClusterChromosomes>>& survivors) {
    survivors = population ;
    survivors.insert(survivors.end(), offsprings.begin(), offsprings.end()) ;  
    sort(survivors.begin(), survivors.end(), [](shared_ptr<ClusterChromosomes>& c1, shared_ptr<ClusterChromosomes>& c2){return c1->fitness<c2->fitness;}) ; 
    for(;survivors.size()>population.size();) survivors.pop_back() ; 

}

void GARouter::Crossover(vector<shared_ptr<ClusterChromosomes>>& population1, 
                        vector<shared_ptr<ClusterChromosomes>>& population2,
                        vector<shared_ptr<ClusterChromosomes>>& offsprings) {
    offsprings.clear() ;

    for(int i=0, p1, p2; offsprings.size()<config.population_size; ++i){
        shared_ptr<ClusterChromosomes> offspring1 = make_shared<ClusterChromosomes>(*population1[i%population1.size()]) ;  
        shared_ptr<ClusterChromosomes> offspring2 = make_shared<ClusterChromosomes>(*population2[i%population2.size()]) ;  
        
        if(crossoverDistrbution(generator)){
            for(int j=0; j<offspring1->size(); ++j){
                p1 = generator() % (1+offspring1->at(j).genotype.size()), p2 = generator() % (1+offspring2->at(j).genotype.size()) ;
                while(p1==p2){ p1 = generator() % (1+offspring1->at(j).genotype.size()) ; p2 = generator() % (1+offspring2->at(j).genotype.size()) ;}
                if(p1>p2) swap(p1, p2) ; 

                for(int k=p1; k<p2; k++) swap(offspring1->at(j).genotype[k], offspring2->at(j).genotype[k]) ;
            }

            offspring1->fitness = -1.0 ;
            offspring2->fitness = -1.0 ;
            if(offsprings.size()<config.population_size) offsprings.push_back(offspring1) ;
            if(offsprings.size()<config.population_size) offsprings.push_back(offspring2) ;
        }else{
            shared_ptr<ClusterChromosomes> offspring1 = population1[i%population1.size()], offspring2 = population2[i%population2.size()] ;
            if(offsprings.size()<config.population_size) offsprings.push_back(offspring1) ;
            if(offsprings.size()<config.population_size) offsprings.push_back(offspring2) ;
        }

    }
}

void GARouter::Crossover(vector<ClusterChromosomes>& population1, vector<ClusterChromosomes>& population2, vector<ClusterChromosomes>& offsprings) {
    offsprings.clear() ;

    for(int i=0, p1, p2; offsprings.size()<config.population_size; ++i){
        ClusterChromosomes offspring1 = population1[i%population1.size()], offspring2 = population2[i%population2.size()] ;
        
        if(crossoverDistrbution(generator)){
            for(int j=0; j<offspring1.size(); ++j){
                p1 = generator() % (1+offspring1[j].genotype.size()), p2 = generator() % (1+offspring2[j].genotype.size()) ;
                while(p1==p2){ p1 = generator() % (1+offspring1[j].genotype.size()) ; p2 = generator() % (1+offspring2[j].genotype.size()) ;}
                if(p1>p2) swap(p1, p2) ; 

                for(int k=p1; k<p2; k++) swap(offspring1[j].genotype[k], offspring2[j].genotype[k]) ;
            }
        }

        if(offsprings.size()<config.population_size) offsprings.push_back(offspring1) ;
        if(offsprings.size()<config.population_size) offsprings.push_back(offspring2) ;
    }
}


void GARouter::Mutation(vector<shared_ptr<ClusterChromosomes>>& population, 
                        vector<shared_ptr<ClusterChromosomes>>& offsprings) {
    offsprings = population ; 
    for(int i=0; i<offsprings.size(); ++i){
        for(int j=0, firstMutate=0; j<offsprings[i]->size(); ++j){
            for(int k=offsprings[i]->at(j).genotype.startIndex; k<offsprings[i]->at(j).genotype.endIndex; ++k){
                if(mutationDistrbution(generator)){ 
                    if(!firstMutate){
                        offsprings[i] = make_shared<ClusterChromosomes>(*offsprings[i]) ; offsprings[i]->fitness = -1.0 ;
                        firstMutate = 1 ; 
                    }
                    offsprings[i]->at(j).genotype[k] = 1-offsprings[i]->at(j).genotype[k] ;
                }
            }
        }
    }
}
void GARouter::Mutation(vector<ClusterChromosomes>& population, vector<ClusterChromosomes>& offsprings){
    offsprings = population ; 

    for(int i=0; i<offsprings.size(); ++i){
        for(int j=0; j<offsprings[i].size(); ++j){
            for(int k=0; k<offsprings[i][j].genotype.size(); ++k){
                if(mutationDistrbution(generator) && offsprings[i][j].genotype[k]!=2) offsprings[i][j].genotype[k] = 1-offsprings[i][j].genotype[k] ;
            }
        }
    }
}

void GARouter::SelectParents(   vector<shared_ptr<ClusterChromosomes>>& population, 
                                vector<shared_ptr<ClusterChromosomes>>& parents1, 
                                vector<shared_ptr<ClusterChromosomes>>& parents2){
    parents1.clear() ; parents2.clear() ; 
    for(int i=0; i<config.population_size; i++){
        vector<shared_ptr<ClusterChromosomes>> candidates ;
        vector<shared_ptr<ClusterChromosomes>>& parents = (i%2) ? parents2 : parents1 ; 
        for(int j=0; j<config.tournament_size; j++) candidates.push_back(population[generator()%population.size()]) ; 

        // 絕對型競爭法
        parents.push_back(*min_element(candidates.begin(), candidates.end(), 
            [](shared_ptr<ClusterChromosomes>& c1, shared_ptr<ClusterChromosomes>& c2){return c1->fitness<c2->fitness;}) ) ; 
    }
}

void GARouter::SelectParents(vector<ClusterChromosomes>& population, vector<ClusterChromosomes>& parents1, vector<ClusterChromosomes>& parents2){
    parents1.clear() ; parents2.clear() ; 
    for(int i=0; i<config.population_size; i++){
        vector<ClusterChromosomes> candidates ;
        vector<ClusterChromosomes>& parents = (i%2) ? parents2 : parents1 ; 
        for(int j=0; j<config.tournament_size; j++) candidates.push_back(population[generator()%population.size()]) ; 

        // 絕對型競爭法
        parents.push_back(*min_element(candidates.begin(), candidates.end(), [](ClusterChromosomes& c1, ClusterChromosomes& c2){return c1.fitness<c2.fitness;}) ) ; 
    }
}

void GARouter::Initial( const vector<tuple<int,int, int, int, int>>& phenotypeInformations, const vector<tuple<int,int,int>>& genotypeInformations, 
                        const vector<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>>& connections,  
                        const vector<ChannelType>& channels, vector<shared_ptr<ClusterChromosomes>>& population){

    population = vector<shared_ptr<ClusterChromosomes>>(config.population_size) ; 
    crossoverDistrbution = bernoulli_distribution(config.cross_prob) ;  
    mutationDistrbution = bernoulli_distribution(config.mut_prob) ;  
    
    for(int i=0; i<population.size(); ++i){
        population[i] = make_shared<ClusterChromosomes>() ;
        ClusterChromosomes& cluster = *population[i] ;
       
        for(int j=0; j<phenotypeInformations.size(); ++j){
            ChromosomeType chromosome = {
                connections[j].first, 
                connections[j].second,
            {
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

            // double r = 5*(j%4) ; 
            // cout << r << "\n" ;
            // for(int k=0; k<chromosome.phenotype.size()-1; ++k){
            //     debugNets.push_back(Net("DUMMY", {{chromosome.phenotype[k]->x+r , chromosome.phenotype[k]->y+r,
            //          chromosome.phenotype[k+1]->x+r , chromosome.phenotype[k+1]->y+r}})) ;
            // }
            // break;
        }
    }
}


void GARouter::Initial( const vector<tuple<int,int, int, int, int>>& phenotypeInformations, const vector<tuple<int,int,int>>& genotypeInformations, 
                        const vector<pair<shared_ptr<ViaNode2>, shared_ptr<ViaNode2>>>& connections,  
                        const vector<ChannelType>& channels, vector<ClusterChromosomes>& population){

    population = vector<ClusterChromosomes>(config.population_size) ; 
    crossoverDistrbution = bernoulli_distribution(config.cross_prob) ;  
    mutationDistrbution = bernoulli_distribution(config.mut_prob) ;  
    
    for(int i=0; i<population.size(); ++i){
        ClusterChromosomes& cluster = population[i] ;
       
        for(int j=0; j<phenotypeInformations.size(); ++j){
            ChromosomeType chromosome = {
                connections[j].first, 
                connections[j].second,
            {
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

            // double r = 5*(j%4) ; 
            // cout << r << "\n" ;
            // for(int k=0; k<chromosome.phenotype.size()-1; ++k){
            //     debugNets.push_back(Net("DUMMY", {{chromosome.phenotype[k]->x+r , chromosome.phenotype[k]->y+r,
            //          chromosome.phenotype[k+1]->x+r , chromosome.phenotype[k+1]->y+r}})) ;
            // }
            // break;
        }
    }
}

bool GenotypeType::operator<(const GenotypeType& genotype){ // A>B
    if(size()!=genotype.size()) throw runtime_error("Unmatched genotype size.") ;
    
    for(int i=startIndex; i<endIndex; ++i){
        if(at(i)==0 && genotype[i]==1) return true ; 
    }
    return false ; 
}

bool GenotypeType::operator>(const GenotypeType& genotype){ // A>B
    if(size()!=genotype.size()) throw runtime_error("Unmatched genotype size.") ;
    
    for(int i=0; i<genotype.size(); ++i){
        if(at(i) && genotype[i]){
            if(at(i)>genotype[i]) return true ; 
        }
    }

    return false ; 
}

double GARouter::CacluteWireLength(ClusterChromosomes& chromosomes){
    double wirtLength = 0.0 ; 
    int p = 2 ;
    for(int i=0; i<chromosomes.size(); ++i){
        for(int j=0; j<chromosomes[i].genotype.size(); ++j){
            if(chromosomes[i].genotype[j]!=p){
                if(!(chromosomes[i].genotype[j]==2 || p==2)) wirtLength += 100 ;
                p = chromosomes[i].genotype[j] ; 
            }
        }
        // for(int j=1; j<chromosomes[i].phenotype.size(); ++j){
        //     auto& tileNode = chromosomes[i].phenotype[j] ; 
        //     auto& prevTileNode = chromosomes[i].phenotype[j-1] ; 
        //     wirtLength += sqrt(pow(prevTileNode->x-tileNode->x, 2) + pow(prevTileNode->y-tileNode->y, 2)) ; 
        // }
    }
    return wirtLength ; 
}
void GARouter::DetermineOrder(ClusterChromosomes& chromosomes, vector<int>& topDownOrders){
    topDownOrders.clear() ; 
    vector<vector<int>> outEdges(chromosomes.size()) ; 
    vector<int> inEdgesCount(chromosomes.size(), 0) ; 
    vector<int> stack1, stack2 ; 

    map<int, vector<int>> groupIndexes ;

    set<int> usedChannels ;
    vector<int> channelOrders ; 

// for(int i=0; i<chromosomes.size(); ++i) topDownOrders.push_back(i) ; 
// return ;
    int conflictCount = chromosomes.size() ; 
    for(int i=0, j=0; i<chromosomes.size(); ){
        while(j<chromosomes.size() && chromosomes[i].phenotype.channelNum==chromosomes[j].phenotype.channelNum) ++j ; 
        for(int k=i; k<j; ++k){
            for(int l=i; l<j; ++l){
                if(k==l) continue ;
                GenotypeType expandedGenotype1=chromosomes[l].genotype, expandedGenotype2=chromosomes[k].genotype ; 
                
                expandedGenotype1.resize(2*chromosomes[l].genotype.size()-1) ; 
                expandedGenotype2.resize(2*chromosomes[k].genotype.size()-1) ; 

                for(int i=0; i<chromosomes[l].genotype.size(); ++i){
                    expandedGenotype1[i*2] = chromosomes[l].genotype[i] ; 
                    if(2*i+1<expandedGenotype1.size()) expandedGenotype1[i*2+1] = chromosomes[l].genotype[i] ; 
                }

                for(int i=0; i<chromosomes[k].genotype.size(); ++i){
                    expandedGenotype2[i*2] = chromosomes[k].genotype[i] ; 
                    if(2*i+1<expandedGenotype2.size()) expandedGenotype2[i*2+1] = chromosomes[k].genotype[i] ; 
                }

                
                for(int i=0; i<expandedGenotype1.size(); ++i){
                    if(expandedGenotype1[i]==1) expandedGenotype1[i] = -2 ;
                    else if(expandedGenotype1[i]==2) expandedGenotype1[i] = 0 ; 
                    else if(expandedGenotype1[i]==0) expandedGenotype1[i] = 2 ; 
                }
                for(int i=0; i<expandedGenotype2.size(); ++i){
                    if(expandedGenotype2[i]==1) expandedGenotype2[i] = -2 ;
                    else if(expandedGenotype2[i]==2) expandedGenotype2[i] = 0 ; 
                    else if(expandedGenotype2[i]==0) expandedGenotype2[i] = 2 ; 
                }
                
                if(chromosomes[l].genotype.startIndex!=0) expandedGenotype1[2*chromosomes[l].genotype.startIndex-1] = (expandedGenotype1[2*chromosomes[l].genotype.startIndex]==2) ? 1 : -1 ;
                if(chromosomes[k].genotype.startIndex!=0) expandedGenotype2[2*chromosomes[k].genotype.startIndex-1] = (expandedGenotype2[2*chromosomes[k].genotype.startIndex]==2) ? 1 : -1 ;

                if(chromosomes[l].genotype.endIndex!=chromosomes[l].genotype.codeLength) expandedGenotype1[2*(chromosomes[l].genotype.endIndex-1)+1] = (expandedGenotype1[2*chromosomes[l].genotype.endIndex-1]==2) ? 1 : -1 ;
                if(chromosomes[k].genotype.endIndex!=chromosomes[k].genotype.codeLength) expandedGenotype2[2*(chromosomes[k].genotype.endIndex-1)+1] = (expandedGenotype2[2*chromosomes[k].genotype.endIndex-1]==2) ? 1 : -1 ;
                
                if( expandedGenotype1>expandedGenotype2 ){
                    outEdges[l].push_back(k) ; // chromosomes[l]比chromosomes[k]上面
                    ++ inEdgesCount[k] ;
                }
            }
        }
        i = j ; 
    }

    
    for(int i=0; i<inEdgesCount.size(); ++i){
        if(!inEdgesCount[i]){
            stack1.push_back(i) ;
            groupIndexes[chromosomes[i].phenotype.channelNum].push_back(i) ; 
        }
    }


    conflictCount -= stack1.size() ;

    for(int i=0; stack1.size() || stack2.size(); ++i){
        vector<int>& currentStack = i%2 ? stack2 : stack1 ;
        vector<int>& nextStack = i%2 ? stack1 : stack2 ;

        for(auto& s : currentStack){
            for(auto& t : outEdges[s]){
                --inEdgesCount[t] ;
                if(!inEdgesCount[t]){
                    nextStack.push_back(t) ;
                    groupIndexes[chromosomes[t].phenotype.channelNum].push_back(t) ; 
                }
            }
        }
        conflictCount -= nextStack.size() ; 
        currentStack.clear() ; 
    }
    
    if(conflictCount!=0) throw runtime_error("Conflict global route!!!") ;

    for(int i=0; i<chromosomes.size(); ++i){
        if(usedChannels.find(chromosomes[i].phenotype.channelNum)==usedChannels.end()){
            usedChannels.insert(chromosomes[i].phenotype.channelNum) ;
            channelOrders.push_back(chromosomes[i].phenotype.channelNum) ;
        }
    }

    sort(channelOrders.begin(), channelOrders.end(), [&chromosomes](int& pa, int& pb)
                    {return chromosomes[0].channelsPtr->at(pa).upChannel[0]->y < chromosomes[0].channelsPtr->at(pb).upChannel[0]->y;}) ;

                
    for(int i=0; i<channelOrders.size(); ++i){
        topDownOrders.insert(topDownOrders.end(), groupIndexes[channelOrders[i]].begin(), groupIndexes[channelOrders[i]].end()) ;
    }
}


double GARouter::CacluteConflictCount(ClusterChromosomes& chromosomes){

    vector<vector<int>> outEdges(chromosomes.size()) ; 
    vector<int> inEdgesCount(chromosomes.size(), 0) ; 
    vector<int> stack1, stack2 ; 

    int conflictCount = chromosomes.size() ; 

    for(int i=0, j=0; i<chromosomes.size(); ){
        while(j<chromosomes.size() && chromosomes[i].phenotype.channelNum==chromosomes[j].phenotype.channelNum) ++j ; 
        for(int k=i; k<j; ++k){
            for(int l=i; l<j; ++l){
                if(k==l) continue ;
                GenotypeType expandedGenotype1=chromosomes[l].genotype, expandedGenotype2=chromosomes[k].genotype ; 
                
                expandedGenotype1.resize(2*chromosomes[l].genotype.size()-1) ; 
                expandedGenotype2.resize(2*chromosomes[k].genotype.size()-1) ; 

                for(int i=0; i<chromosomes[l].genotype.size(); ++i){
                    expandedGenotype1[i*2] = chromosomes[l].genotype[i] ; 
                    if(2*i+1<expandedGenotype1.size()) expandedGenotype1[i*2+1] = chromosomes[l].genotype[i] ; 
                }

                for(int i=0; i<chromosomes[k].genotype.size(); ++i){
                    expandedGenotype2[i*2] = chromosomes[k].genotype[i] ; 
                    if(2*i+1<expandedGenotype2.size()) expandedGenotype2[i*2+1] = chromosomes[k].genotype[i] ; 
                }

                
                for(int i=0; i<expandedGenotype1.size(); ++i){
                    if(expandedGenotype1[i]==1) expandedGenotype1[i] = -2 ;
                    else if(expandedGenotype1[i]==2) expandedGenotype1[i] = 0 ; 
                    else if(expandedGenotype1[i]==0) expandedGenotype1[i] = 2 ; 
                }
                for(int i=0; i<expandedGenotype2.size(); ++i){
                    if(expandedGenotype2[i]==1) expandedGenotype2[i] = -2 ;
                    else if(expandedGenotype2[i]==2) expandedGenotype2[i] = 0 ; 
                    else if(expandedGenotype2[i]==0) expandedGenotype2[i] = 2 ; 
                }
                
                if(chromosomes[l].genotype.startIndex!=0) expandedGenotype1[2*chromosomes[l].genotype.startIndex-1] = (expandedGenotype1[2*chromosomes[l].genotype.startIndex]==2) ? 1 : -1 ;
                if(chromosomes[k].genotype.startIndex!=0) expandedGenotype2[2*chromosomes[k].genotype.startIndex-1] = (expandedGenotype2[2*chromosomes[k].genotype.startIndex]==2) ? 1 : -1 ;

                if(chromosomes[l].genotype.endIndex!=chromosomes[l].genotype.codeLength) expandedGenotype1[2*(chromosomes[l].genotype.endIndex-1)+1] = (expandedGenotype1[2*chromosomes[l].genotype.endIndex-1]==2) ? 1 : -1 ;
                if(chromosomes[k].genotype.endIndex!=chromosomes[k].genotype.codeLength) expandedGenotype2[2*(chromosomes[k].genotype.endIndex-1)+1] = (expandedGenotype2[2*chromosomes[k].genotype.endIndex-1]==2) ? 1 : -1 ;
                
                if(showExpand){
                    cout << chromosomes[l].genotype << "|\n" ;
                    cout << expandedGenotype1 << "##\n" ;
                    cout << chromosomes[k].genotype << "|\n" ;
                    cout << expandedGenotype2 << "##\n" ;
                    cout << (expandedGenotype1>expandedGenotype2) << "\n" ;
                }

                if( expandedGenotype1>expandedGenotype2 ){
                    outEdges[l].push_back(k) ; // chromosomes[l]比chromosomes[k]上面
                    ++ inEdgesCount[k] ;
                }
            }
        }
        i = j ; 
    }

    for(int i=0; i<inEdgesCount.size(); ++i) if(!inEdgesCount[i]) stack1.push_back(i) ;
    conflictCount -= stack1.size() ;

    for(int i=0; stack1.size() || stack2.size(); ++i){
        vector<int>& currentStack = i%2 ? stack2 : stack1 ;
        vector<int>& nextStack = i%2 ? stack1 : stack2 ;

        for(auto& s : currentStack){
            for(auto& t : outEdges[s]){
                --inEdgesCount[t] ;
                if(!inEdgesCount[t]) nextStack.push_back(t) ;
            }
        }
        conflictCount -= nextStack.size() ; 
        currentStack.clear() ; 
    }

    return conflictCount ; 
}

double GARouter::CacluteCapacityValue(ClusterChromosomes& chromosomes){
    unordered_map<shared_ptr<int>, int> originalCapacity ;

    double incompatibleCount = 0.0 ; 
    for(int i=0; i<chromosomes.size(); ++i){
        for(int j=0; j<chromosomes[i].phenotype.size(); ++j){
            if(originalCapacity.find(chromosomes[i].phenotype[j].capacity)==originalCapacity.end()) originalCapacity[chromosomes[i].phenotype[j].capacity] = *chromosomes[i].phenotype[j].capacity ;
            -- (*chromosomes[i].phenotype[j].capacity) ;
        }
    }

    for(auto& [capacity, oldCapacity] : originalCapacity){
        if((*capacity)<0) incompatibleCount -= *capacity ;
        *(capacity) = oldCapacity ; 
    }
    return incompatibleCount ; 
//     double incompatibleCount = 0.0 ; 
//     for(int i=0; i<chromosomes.size(); ++i){
//         for(int j=0; j<chromosomes[i].phenotype.size(); ++j){
//             ++ (*chromosomes[i].phenotype[j].currentCapacity) ;
//         }
//     }

//    for(int i=0; i<chromosomes.size(); ++i){
//         for(int j=0; j<chromosomes[i].phenotype.size(); ++j){
//             if(*chromosomes[i].phenotype[j].currentCapacity < 0) incompatibleCount -= *chromosomes[i].phenotype[j].currentCapacity ;
//             *chromosomes[i].phenotype[j].currentCapacity = 0 ;
//         }
//     }
//     return incompatibleCount ; 
}
double GARouter::Fitness(ClusterChromosomes& chromosomes){
    double length = CacluteWireLength(chromosomes) ;
    double incapaciableCount = CacluteCapacityValue(chromosomes) ;
    double conflictCount = CacluteConflictCount(chromosomes) ; 
    return config.alpha*length + config.beta*incapaciableCount + config.gamma*conflictCount ; 
}

void GARouter::Evaluate(vector<shared_ptr<ClusterChromosomes>>& population){
    for(auto& cluster : population){
        if(cluster->fitness<0) cluster->fitness = Fitness(*cluster) ; 
    }
}
void GARouter::Evaluate(vector<ClusterChromosomes>& population){
    for(auto& cluster : population){
        cluster.fitness = Fitness(cluster) ; 
    }
}


void GARouter::EvaluateAll(vector<shared_ptr<ClusterChromosomes>>& population,bool showAverage){
    Evaluate(population) ; 
    // double sum = 0.0 ; 
    // for(auto& c : population) sum += c.fitness ;

    // cout << sum/population.size() << "\n" ;
}


void GARouter::EvaluateAll(vector<ClusterChromosomes>& population,bool showAverage){
    Evaluate(population) ; 
    // double sum = 0.0 ; 
    // for(auto& c : population) sum += c.fitness ;

    // cout << sum/population.size() << "\n" ;
}
