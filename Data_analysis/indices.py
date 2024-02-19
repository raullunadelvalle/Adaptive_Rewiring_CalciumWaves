# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 10:17:02 2020

@author: JiaLi
"""
import numpy as np
import networkx as nx
from scripts import basic_comp as bcomp

# from adjactent matrix to graph
def convert_from_adj2networkX(A,weighted=True):

    edges_ind = np.where(A>0)
    num_edges = len(edges_ind[0])

    G = nx.DiGraph()   # DiGraph
    G.add_nodes_from(np.arange(A.shape[0]))

    edges_list = list()
    if weighted:
        for ind in np.arange(num_edges):
            edge_pair_w = (edges_ind[1][ind],edges_ind[0][ind],1/A[edges_ind[0][ind],edges_ind[1][ind]]) #distance = 1/weight
            edges_list.append(edge_pair_w)    
        
        G.add_weighted_edges_from(edges_list)
    else:
        for ind in np.arange(num_edges):
            edge_pair = (edges_ind[1][ind],edges_ind[0][ind])
            edges_list.append(edge_pair)
        
        G.add_edges_from(edges_list)    
        
        
    return G

# compute path_length and connectivity
# efficiency: newman, 2003, The Structure and Function of Complex Networks
def path_length_connectedness(A):
    G =  convert_from_adj2networkX(A,True)
    len_paths = dict(nx.all_pairs_dijkstra_path_length(G))
    connected_pairs = np.zeros(A.shape)
    connected_pairs_length = np.zeros(A.shape)
    # path from s to t, s differs from t
    for s in len_paths.keys():
        for t in len_paths[s].keys():
            if(len_paths[s][t]>0):
                connected_pairs[s][t] = 1
                connected_pairs_length[s][t] = len_paths[s][t]
    
    # path from s to s
    for i in range(A.shape[0]):
        self_loop = np.where((connected_pairs[i,:]+connected_pairs[:,i])>1)[0]
        if(len(self_loop)>0):
            connected_pairs[i,i] = 1
            connected_pairs_length[i,i] = np.min((connected_pairs_length[i,:]+connected_pairs_length[:,i])[self_loop])
        
    connected_pairs_length[connected_pairs_length==0] = np.inf

    return  connected_pairs, connected_pairs_length

# number of hubs
def hub_number(A,thresh,binary_flag,axisUsed=0):
    if binary_flag==True:
        deg_u = np.sum(A>0,axis=axisUsed)
        deg_n = np.sum(A>0,axis=1-axisUsed)
    else:
        deg_u = np.sum(A,axis=axisUsed)
        deg_n = np.sum(A,axis=1-axisUsed)
        
    num_hubs = len(np.where((deg_u>=thresh)&(deg_n>0))[0])
    return num_hubs

# convergent-divergent unit
def cd_pairs(A, connected_pairs, thresh):
    
    Topo_A = 1.0*(A>0)
    deg_in = np.sum(Topo_A,1)
    deg_out = np.sum(Topo_A,0)    
    candCon = np.where((deg_in>=thresh)&(deg_out>0))[0]
    candDiv = np.where((deg_out>=thresh)&(deg_in>0))[0]
    connected_pairs_noloop = connected_pairs.copy()
    np.fill_diagonal(connected_pairs_noloop,0)
    
    if (len(candCon)>0)&(len(candDiv)>0):
        cd_pairs_temp = np.where(connected_pairs_noloop[np.ix_(candCon,candDiv)]>0)
        cHubs = candCon[cd_pairs_temp[0]]
        dHubs = candDiv[cd_pairs_temp[1]]
        connected_cd_pairs = np.stack((cHubs,dHubs), axis=0)
    else:
        connected_cd_pairs = np.array([[],[]])
        
    return connected_cd_pairs

# source and target nodes
def sour_targ_number(connected_cd_pairs, connected_pairs):
    cHubs = connected_cd_pairs[0]
    dHubs = connected_cd_pairs[1]

    connected_pairs_noloop = connected_pairs.copy()
    np.fill_diagonal(connected_pairs_noloop,0)
    targ_nodes = connected_pairs_noloop[dHubs,:]
    sour_nodes = connected_pairs_noloop[:,cHubs]
    for i in range(len(dHubs)):
        targ_nodes[i,cHubs[i]] = 0
        sour_nodes[dHubs[i],i] = 0
        
    num_targ = np.sum(targ_nodes,axis=1)
    num_sour = np.sum(sour_nodes,axis=0)
        
    total_nodes = 1.0*(targ_nodes+np.transpose(sour_nodes)>0)
    overlap_nodes = 1.0*(targ_nodes+np.transpose(sour_nodes)>1)
    
    num_total = np.sum(total_nodes,axis=1)
    num_overlap = np.sum(overlap_nodes,axis=1)
    prop_overlap = num_overlap/num_total
    
    return num_targ, num_sour, prop_overlap

# size and density of intermediate subgraphs
def intermediate_subgraphs(A,connected_cd_pairs):
    Topo_A = 1.0*(A>0)
    deg_in = np.sum(Topo_A,1)
    deg_out = np.sum(Topo_A,0)
    num_cd = connected_cd_pairs.shape[1]
    n = A.shape[0]
    
    interm_size = []
    interm_density = []
    periph_density = []
    
    for i in range(num_cd):
        c = connected_cd_pairs[0,i]
        d = connected_cd_pairs[1,i]
        interm = []
        
        B = Topo_A.copy()
        B[c,:] = 0 #cut in-links of c
        B[:,d] = 0 #cut out-links of d
        G = convert_from_adj2networkX(B,True)           
        
        nodes = set(np.arange(n)).difference(set([c,d]))        #all nodes except c and d
        for k in nodes:
            if nx.has_path(G,c,k) & nx.has_path(G,k,d): # c->k->d exist
                interm.append(k)

        interm_size.append(len(interm))
        periph_nodes = list(nodes.difference(set(interm)))
        if (len(interm)>1):
            interm_density.append(np.sum(B[np.ix_(interm, interm)])/(len(interm)*(len(interm)-1)))   
        else:#RL adds
            interm_density.append(np.nan)#RL adds
         
        if (len(periph_nodes)>1):
            periph_density.append(np.sum(B[np.ix_(periph_nodes, periph_nodes)])/(len(periph_nodes)*(len(periph_nodes)-1)))
        else:#RL adds
            periph_density.append(np.nan)#RL adds
                            
    return interm_size, interm_density,periph_density

#A_inv_paths[i,j] is the inverse path from i to j
def get_inv_path_all(G):
    
    
    len_paths = dict(nx.all_pairs_dijkstra_path_length(G))
        
    nodes = len(G.nodes)
    A_inv_paths = np.zeros((nodes,nodes))
    
    for n in range(nodes): 
        A_inv_paths[n,n] = np.nan #the inverted distance of a node to itself is infinite
    
    
    for tN in range(nodes): #take each node to test
        for k in len_paths[tN].keys():
            
            if k != tN:            
                A_inv_paths[tN,k] = 1/len_paths[tN][k]         
        
    return A_inv_paths

#gives you back the dictionary connections where at each entry (same keys as the combinations of parameters)
# you have an array that has the number of hubs, for different defs of hubs 
#(starting from min_con and incrementing by one until you reach max_con)


# modularity
# Leicht, Elizabeth A., and Mark EJ Newman. 2008, "Community structure in directed networks."

def modularity(A_matrices):
    Q = {}
    communitiesDict = {}
    for key in A_matrices.keys():
        Q[key] = {}
        for k in A_matrices[key].keys():
            Q[key][k], communitiesDict[key][k] = directed_modularity.getModularityIndex(A_matrices[key][k])

    return Q

# In[clustering]
# Fagiolo, 2007, Clustering in complex directed networks

def weighted_directed_clustering(A):
    G = nx.from_numpy_matrix(A.T,create_using=nx.DiGraph)
    clustering = nx.average_clustering(G, weight='weight')
    
    return clustering

# In[small world]
# small world indices: C, E, S

def small_world_indices(A_matrices):
    clustering = {}
    efficiency = {}
    small_world = {}

    for key in A_matrices.keys():
        clustering[key] = {}
        efficiency[key] = {}
        small_world[key] = {}
        for k in A_matrices[key].keys():
            clustering[key][k] = weighted_directed_clustering(A_matrices[key][k])
            efficiency[key][k] = weighted_directed_efficiency(A_matrices[key][k])
            small_world[key][k] = clustering[key][k]*efficiency[key][k]
            
    return small_world, clustering, efficiency

def small_world_indices0(A_matrices):
    clustering = {}
    efficiency = {}
    small_world = {}

    for key in A_matrices.keys():
        clustering[key] = weighted_directed_clustering(A_matrices[key])
        efficiency[key] = weighted_directed_efficiency(A_matrices[key])
        small_world[key] = clustering[key]*efficiency[key]
            
    return small_world, clustering, efficiency


# In[null model]

# generate random adjacent matrix randomizing topology or weights

def generate_match_adj(A, flag_rnd, **kwargs):
    A_rand = 1.0*(A>0)
    
    if (flag_rnd == 'topology'):
        for (key, value,) in kwargs.items():
            if key == "swaps":
                swaps = value
            else:
                print('Please specify the number of swaps')
                return -1       
        
        n_edges = np.sum(A>0)
        for i in range(swaps):
            non_zero_ind = np.where(A_rand>0)
            common_nodes = [A.shape[0]+1]
            flag_1 = 999
            flag_2 = 999
            while (len(common_nodes)!=0 or (flag_1+flag_2)>0):
                swap_edges  = np.random.choice(n_edges,size=2,replace=False)
                edge_1 = (non_zero_ind[0][swap_edges[0]],non_zero_ind[1][swap_edges[0]])
                edge_2 = (non_zero_ind[0][swap_edges[1]],non_zero_ind[1][swap_edges[1]])
                common_nodes = list(set(edge_1).intersection(edge_2))
                flag_1 = A_rand[edge_1[0],edge_2[1]]
                flag_2 = A_rand[edge_2[0],edge_1[1]]
            A_rand[edge_1[0],edge_2[1]] = 1
            A_rand[edge_2[0],edge_1[1]] =1
            A_rand[edge_1] = 0
            A_rand[edge_2] = 0
            
    if (flag_rnd == 'weights'):
        weight_list = A[np.where(A>0)]
        np.random.shuffle(weight_list)
        A_rand[np.where(A_rand>0)] = weight_list
        
    return A_rand

# In[standardize]
#Crand & Lrand

def random_clustering_and_efficiency(A_rand):
    m = len(A_rand)
    rand_clustering = np.zeros(m)
    rand_efficiency = np.zeros(m)

    for k in range(m):
        rand_clustering[k] = weighted_directed_clustering(A_rand[k])
        rand_efficiency[k] = weighted_directed_efficiency(A_rand[k])
    Crand = np.mean(rand_clustering)
    Erand = np.mean(rand_efficiency)
    
    return Crand, Erand

# In[motifs]
#<script src="https://gist.github.com/tpoisot/8582648.js"></script>

def motif_counter(Topo_A, n_motif, motifs):

    mcount = np.repeat(0,len(motifs))

    G = nx.from_numpy_matrix(Topo_A.T, create_using=nx.DiGraph)
    nodes = G.nodes()
    
    tuples_set = itertools.combinations(nodes,n_motif)

    for tuples in list(tuples_set):
        #print(tuples)
        sub_gr = G.subgraph(tuples)
        for ind, key in enumerate(motifs):
            mot_match = nx.is_isomorphic(sub_gr, motifs[key])
            if mot_match:
                mcount[ind] += 1

    return mcount


def topo_null_model(Topo_A, swaps):
    n_edges = np.sum(Topo_A>0)
    A_rand = Topo_A.copy()
    for i in range(swaps):
        non_zero_ind = np.where(A_rand>0)
        common_nodes = [A_rand.shape[0]+1]
        flag_1 = 999
        flag_2 = 999
        while (len(common_nodes)!=0 or (flag_1+flag_2)>0):
            swap_edges  = np.random.choice(n_edges,size=2,replace=False)
            edge_1 = (non_zero_ind[0][swap_edges[0]],non_zero_ind[1][swap_edges[0]])
            edge_2 = (non_zero_ind[0][swap_edges[1]],non_zero_ind[1][swap_edges[1]])
            common_nodes = list(set(edge_1).intersection(edge_2))
            flag_1 = A_rand[edge_1[0],edge_2[1]]
            flag_2 = A_rand[edge_2[0],edge_1[1]]
        A_rand[edge_1[0],edge_2[1]] = 1
        A_rand[edge_2[0],edge_1[1]] =1
        A_rand[edge_1] = 0
        A_rand[edge_2] = 0
        
    return A_rand