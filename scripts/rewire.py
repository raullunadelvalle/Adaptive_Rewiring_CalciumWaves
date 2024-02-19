# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 12:13:46 2020

@author: JiaLi
"""

import numpy as np
import sys
sys.path.append("..")
from scripts import basic_comp as bcomp
from scripts import initialization_params

import seaborn as sns
import matplotlib.pylab as plt

import networkx as nx


# initialize a random directed weighted network
def initial_directed_network(n_nodes, edges, weight_distribution,**kwargs):
    rand_adj_matrix = bcomp.generate_rand_adj(n_nodes, edges, weight_distribution, **kwargs)
    
    nodes_coord = bcomp.generate_nodes_coord(n_nodes)
    
    return nodes_coord, rand_adj_matrix

# replicate Rentzeperis et al., 2022
def rewRand(rand_adj_matrix, tau, p, p_in, rewirings):
    
    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    
    for iter_i in range(rewirings):
        r_in = np.random.random_sample()
        if r_in<p_in:
            deg_in = np.sum(A > 0, axis=1, keepdims=False,)
            nodes_receiving = np.where((deg_in > 0) & (deg_in < n- 1))
            if len(nodes_receiving)==0:
                print("All nodes have either 0 or (n-1) in-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_receiving[0], n, flag='in')

            flag_alg = 'consensus'
            flag = 'in'
            
            r = np.random.random_sample()    
            j_add, j_minus = bcomp.rewire_random(r,p,tau,A,v,U_cp_v,U_nc_v,flag,flag_alg)
            
            A[v, j_add,] = A[v, j_minus,]
            A[v, j_minus, ] = 0
        
        else:
            deg_out = np.sum(A > 0, axis=0, keepdims=False,)
            nodes_sending = np.where((deg_out > 0) & (deg_out < n - 1))
            if len(nodes_sending)==0:
                print("All nodes have either 0 or (n-1) out-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_sending[0], n, flag='out')

            flag_alg = 'advection'
            flag = 'out'
            
            r = np.random.random_sample()    
            j_add, j_minus = bcomp.rewire_random(r,p,tau,A,v,U_cp_v,U_nc_v,flag,flag_alg)
        
            A[j_add, v, ] = A[j_minus, v, ]
            A[j_minus, v, ] = 0
            
    return A

# substitute distance-based rewiring for random rewiring
def rewDist(vertices_coord, rand_adj_matrix, tau, p, p_in, rewirings):
    
    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    D = bcomp.compute_distance_matrix(vertices_coord)
    
    for iter_i in range(rewirings):
        r_in = np.random.random_sample()
        if r_in<p_in:
            deg_in = np.sum(A > 0, axis=1, keepdims=False,)
            nodes_receiving = np.where((deg_in > 0) & (deg_in < n- 1))
            if len(nodes_receiving)==0:
                print("All nodes have either 0 or (n-1) in-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_receiving[0], n, flag='in')

            flag_alg = 'consensus'
            flag = 'in'
            
            r = np.random.random_sample()    
            j_add, j_minus = bcomp.rewire_distance(r,p,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg)
            
            A[v, j_add,] = A[v, j_minus,]
            A[v, j_minus, ] = 0
        
        else:
            deg_out = np.sum(A > 0, axis=0, keepdims=False,)
            nodes_sending = np.where((deg_out > 0) & (deg_out < n - 1))
            if len(nodes_sending)==0:
                print("All nodes have either 0 or (n-1) out-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_sending[0], n, flag='out')

            flag_alg = 'advection'
            flag = 'out'
            
            r = np.random.random_sample()    
            j_add, j_minus = bcomp.rewire_distance(r,p,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg)
        
            A[j_add, v, ] = A[j_minus, v, ]
            A[j_minus, v, ] = 0
            
    return A





# Rewiring that combines diffussion and calcium waves/giant depolarizing potentials (i.e. activity Boosts)
def rewBoost(vertices_coord, rand_adj_matrix, tau, pDist, p_in, rewirings, pBoost,boost_node, booster):

    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    D = bcomp.compute_distance_matrix(vertices_coord)
    
 
    for iter_i in range(rewirings):
 
        r_in = np.random.random_sample()       
        if r_in<p_in:
            deg_in = np.sum(A > 0, axis=1, keepdims=False,)
            nodes_receiving = np.where((deg_in > 0) & (deg_in < n- 1))


            if len(nodes_receiving)==0:
                print("All nodes have either 0 or (n-1) in-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_receiving[0], n, flag='in')

            flag_alg = 'consensus'
            flag = 'in'
            
            r = np.random.random_sample()    
            j_add, j_minus = bcomp.rewire_distance_boost(r,pDist,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg, pBoost,boost_node, booster)
            
            A[v, j_add,] = A[v, j_minus,]
            A[v, j_minus, ] = 0
        
        else:
            deg_out = np.sum(A > 0, axis=0, keepdims=False,)
            nodes_sending = np.where((deg_out > 0) & (deg_out < n - 1))

              
            if len(nodes_sending)==0:
                print("All nodes have either 0 or (n-1) out-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_sending[0], n, flag='out')

            flag_alg = 'advection'
            flag = 'out'
            
            r = np.random.random_sample()  
            j_add, j_minus = bcomp.rewire_distance_boost(r,pDist,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg, pBoost,boost_node, booster)
        
            A[j_add, v, ] = A[j_minus, v, ]
            A[j_minus, v, ] = 0
               
  
         
    return A






# only rewire out-links
def run_out_dynamics(vertices_coord, rand_adj_matrix, p, q, rewirings, tau, vector_field, flag_alg,):

    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    
    D = bcomp.compute_distance_matrix(vertices_coord)
    W = bcomp.compute_wave_matrix(vertices_coord, vector_field, D, flag='out')
    
    for iter_i in range(rewirings):
        deg_out = np.sum(A > 0, axis=0, keepdims=False,)
        nodes_sending = np.where((deg_out > 0) & (deg_out < n - 1))
        if len(nodes_sending)==0:
            print("All nodes have either 0 or (n-1) out-degree.")
            return A
        
        v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_sending[0], n, flag='out')

        r = np.random.random_sample()       
        j_add, j_minus = bcomp.rewire_out_principles(r,p,q,D,W,tau,A,v,U_cp_v,U_nc_v,flag_alg)
        
        A[j_add, v, ] = A[j_minus, v, ]
        A[j_minus, v, ] = 0

    return A

# only rewire in-links
def run_in_dynamics(vertices_coord, rand_adj_matrix, p, q, rewirings, tau, vector_field, flag_alg,):

    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    
    D = bcomp.compute_distance_matrix(vertices_coord)
    W = bcomp.compute_wave_matrix(vertices_coord, vector_field, D, flag='in')
    
    for iter_i in range(rewirings):
        deg_in = np.sum(A > 0, axis=1, keepdims=False,)
        nodes_receiving = np.where((deg_in > 0) & (deg_in < n- 1))
        if len(nodes_receiving)==0:
            print("All nodes have either 0 or (n-1) in-degree.")
            return A
        
        v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_receiving[0], n, flag='in')

        r = np.random.random_sample()       
        j_add, j_minus = bcomp.rewire_in_principles(r,p,q,D,W,tau,A,v,U_cp_v,U_nc_v,flag_alg)
            
        A[v, j_add,] = A[v, j_minus,]
        A[v, j_minus, ] = 0

    return A
    
# either run advection or consensus at each iteration
def run_dynamics_advection_consensus_sequence(vertices_coord, rand_adj_matrix, p, q, rewirings, tau, vector_field, p_in):
    
    n = rand_adj_matrix.shape[0]
    A = rand_adj_matrix.copy()
    
    D = bcomp.compute_distance_matrix(vertices_coord)
    W_out = bcomp.compute_wave_matrix(vertices_coord, vector_field, D, flag='out')
    W_in = bcomp.compute_wave_matrix(vertices_coord, vector_field, D, flag='in')
    
    for iter_i in range(rewirings):
        r_in = np.random.random_sample()
        if r_in<p_in:
            deg_in = np.sum(A > 0, axis=1, keepdims=False,)
            nodes_receiving = np.where((deg_in > 0) & (deg_in < n- 1))
            if len(nodes_receiving)==0:
                print("All nodes have either 0 or (n-1) in-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_receiving[0], n, flag='in')

            r = np.random.random_sample()    
            flag_alg='consensus'
            j_add, j_minus = bcomp.rewire_in_principles(r,p,q,D,W_in,tau,A,v,U_cp_v,U_nc_v,flag_alg)
            
            A[v, j_add,] = A[v, j_minus,]
            A[v, j_minus, ] = 0
        
        else:
                
            deg_out = np.sum(A > 0, axis=0, keepdims=False,)
            nodes_sending = np.where((deg_out > 0) & (deg_out < n - 1))
            if len(nodes_sending)==0:
                print("All nodes have either 0 or (n-1) out-degree.")
                return A
        
            v, U_cp_v, U_nc_v = bcomp.choose_rewire_vertex(A, nodes_sending[0], n, flag='out')

            r = np.random.random_sample()    
            flag_alg='advection'
            j_add, j_minus = bcomp.rewire_out_principles(r,p,q,D,W_out,tau,A,v,U_cp_v,U_nc_v,flag_alg)
        
            A[j_add, v, ] = A[j_minus, v, ]
            A[j_minus, v, ] = 0
            
    return A
