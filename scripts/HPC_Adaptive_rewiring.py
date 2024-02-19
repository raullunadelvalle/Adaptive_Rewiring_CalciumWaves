import numpy as np
import pickle
import rewire, indices, directed_modularity
from scripts import  initialization_params

import sys
sys.path.append("..") #sys.path.append(".")

f = open('./Output/initials_nor.pckl', 'rb')
vertices_coord, initial_matrix = pickle.load(f)
f.close()
K = len(vertices_coord)

# parameters #################################################
Tau = [1]

p_in_list = [0.5] 
p_dist = list(map(lambda x: x/10,range(11))) #0 distance rewiring never occurs / 1: distance rewiring always occurs whenever there is no boost
p_boost = list(map(lambda x: x/10,range(11))) #calcium wave ocurrence probability

rewirings = 4000
Rep = 1



booster = 1


    
# Rewiring starts here ################################################
A_matrices = {}
boost_node = {}

for k in range(K): #each k is an initial random network
    boost_node[k] = initialization_params.boost_node_selector(initial_matrix[k])


for tau in Tau:
    key_t = 'tau='+str(tau)
    A_matrices[tau] = {}


    for p_in in p_in_list:
        A_matrices[tau][p_in] = {}
        
        
        for pDist in p_dist:
            A_matrices[tau][p_in][pDist] = {}
            for pBoost in p_boost:
              
                
                print(pDist+pBoost)
                if (pDist+pBoost<=1):
                    A_matrices[tau][p_in][pDist][pBoost] = {}

                    
                    for k in range(K): #each k is an initial random network
                        A_matrices[tau][p_in][pDist][pBoost][k] = {}
 
    
                        for r in range(Rep): #each r is a repetition over an initial random network
                            A_matrices[tau][p_in][pDist][pBoost][k][r] = rewire.rewBoost(vertices_coord[k], initial_matrix[k], tau, pDist, p_in, rewirings, pBoost,boost_node[k], booster)



f = open('./Output/simulation_params.pckl', 'wb')
pickle.dump([Tau,p_in_list,p_dist,p_boost,rewirings,K,Rep,booster], f)
f.close()

fname = './Output/boosted_nodes.pckl'
f = open(fname, 'wb')
pickle.dump(boost_node, f)
f.close()

fname = './Output/boost_nor.pckl'
f = open(fname, 'wb')
pickle.dump(A_matrices, f)
f.close()

fname = './Output/boost_nor_MatPlusCoords.pckl'
f = open(fname, 'wb')
pickle.dump([vertices_coord,A_matrices], f) 
f.close()

#####################################################################