# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 16:24:35 2020

@author: JiaLi
"""

import numpy as np
from scipy import linalg
from scipy.spatial.distance import squareform, pdist

import sys
sys.path.append("..")

# generate nodes coordinators
def generate_nodes_coord(n_nodes):
    phi = np.random.random(size=n_nodes) * 2 * np.pi
    r2 = np.random.random(size=n_nodes)
    x = np.sqrt(r2)*np.cos(phi)
    y = np.sqrt(r2)*np.sin(phi)
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    nodes_coord = np.concatenate((x,y),axis=1)
    return nodes_coord

# generate random directed adj
def generate_rand_adj(n_nodes, edges, weight_distribution,**kwargs):

    for (key, value,) in kwargs.items():
        if key == "mu":
            mu = value
        elif key == "sig":
            sig = value

    # set constant
    EPSILON = 0.05

    # set the max number of network edges
    max_connections = int(n_nodes * (n_nodes - 1))
    print("max connections:", max_connections,)
    print("weight distribution:", weight_distribution,)
    
    if edges > max_connections or edges < 0:
        print("Edge number out of range")
        return -1
        
    print("Generating random adjacency matrix ...")

    # sample weights from a lognormal distribution
    if weight_distribution == "lognormal":
        rand_weights = np.random.lognormal(mean=mu, sigma=sig, size=edges,)

    # ... from a normal distribution
    elif weight_distribution == "normal":
        rand_weights = np.random.normal(loc=mu, scale=sig, size=edges,)
        ind = np.where(rand_weights < 0)
        rand_weights[ind] = EPSILON

        # ... from a binary distribution
    elif weight_distribution == "binary":
        rand_weights = np.ones(edges)

    print("Weights generation: random weights: type:{}, ndim: {}, shape: {}, mean: {}, std: {}".format(type(rand_weights),rand_weights.ndim,rand_weights.shape,np.mean(rand_weights),np.std(rand_weights),))

    # Normalize weights such that their sum equals the number of edges
    if (weight_distribution == "normal") | (weight_distribution == "lognormal"):
        norm_factor = edges / np.sum(rand_weights)
        norm_rand_weights = rand_weights * norm_factor
        # [print info on normalized weights for validation]
        print("Weights normalization: normalized random weights: type:{}, ndim: {}, shape: {}, min: {}, max: {},mean: {},std: {},factor: {}".format(type(norm_rand_weights),norm_rand_weights.ndim,norm_rand_weights.shape,norm_rand_weights.min(),norm_rand_weights.max(),np.mean(norm_rand_weights),np.std(norm_rand_weights),norm_factor))
    else:
        norm_rand_weights = rand_weights

    # Get the indices of 1s of a matrix the same size as A with 1s everywhere except in the diagonal
    Aones = np.ones((n_nodes, n_nodes)) - np.eye(n_nodes)
    ind = np.where(Aones)

    # Pick a random sample of those indices (# edges)
    rand_max_con = np.random.permutation(max_connections)
    rand_edges_ind = (ind[0][rand_max_con[:edges]],ind[1][rand_max_con[:edges]],)

    # build the adjacency matrix w/ those indices
    adj_matx = np.zeros((n_nodes, n_nodes))
    adj_matx[rand_edges_ind] = norm_rand_weights
    print("Adjacent matrix: type:{}, ndim: {}, shape: {}, min: {}, max: {}".format(type(adj_matx), adj_matx.ndim, adj_matx.shape, adj_matx.min(), adj_matx.max(), ))

    print("Completed")
    return adj_matx




# consensus kernel
def compute_consensus_kernel(weight_matx, tau):

    # estimate the in degree Laplacian
    Din = np.diag(np.sum(weight_matx, axis=1))
    Lin = Din - weight_matx

    # calculate the consensus kernel
    kernel = linalg.expm(-tau * Lin)

    return kernel

# advection kernel
def compute_advection_kernel(weight_matx, tau):
    
    # estimate the out degree Laplacian
    Dout = np.diag(np.sum(weight_matx, axis=0))
    Lout = Dout - weight_matx

    # calculate the consensus kernel
    kernel = linalg.expm(-tau * Lout)

    return kernel




# matrix of principle 1: distance
def compute_distance_matrix(vertices_coord):
    D = squareform(pdist(vertices_coord))
    return D

# matrix of principle 2: wave
def compute_wave_matrix(vertices_coord, vector_field, D, flag):
    
    denom_cos_theta = D*(linalg.norm(vector_field,axis=1)[:,np.newaxis])
    np.fill_diagonal(denom_cos_theta,1)
        
    if flag == 'out': 
        nom_cos_theta_out = np.sum(-(vertices_coord[:,np.newaxis,:]-vertices_coord)*vector_field[:,np.newaxis],axis=2)
        np.fill_diagonal(nom_cos_theta_out,1)
        W_out = (1-nom_cos_theta_out/denom_cos_theta).T
        return W_out
    elif flag == 'in':
        nom_cos_theta_in = np.sum((vertices_coord[:,np.newaxis,:]-vertices_coord)*vector_field[:,np.newaxis],axis=2)
        np.fill_diagonal(nom_cos_theta_in,1)
        W_in = (1-nom_cos_theta_in/denom_cos_theta)
        return W_in
    else:
        nom_cos_theta_out = np.sum(-(vertices_coord[:,np.newaxis,:]-vertices_coord)*vector_field[:,np.newaxis],axis=2)
        nom_cos_theta_in = np.sum((vertices_coord[:,np.newaxis,:]-vertices_coord)*vector_field[:,np.newaxis],axis=2)
        np.fill_diagonal(nom_cos_theta_in,1)
        np.fill_diagonal(nom_cos_theta_out,1)
        W_out = (1-nom_cos_theta_out/denom_cos_theta).T
        W_in = (1-nom_cos_theta_in/denom_cos_theta)
        return W_in, W_out


# choose the vertex that waited to be rewired for degree, and the sets of nodes that are children/not children/parents/not parents of this vertex 
def choose_rewire_vertex(A, nodes_receiving,n,flag):
    v = np.random.choice(nodes_receiving)
    all_vertices_ind = np.arange(n)
    no_v_vertices_ind = np.delete(all_vertices_ind, v)
    if flag=='in':
        U_cp_v = no_v_vertices_ind[ np.where(A[v, no_v_vertices_ind, ]>0)[0]]
        U_nc_v = no_v_vertices_ind[np.where(A[v, no_v_vertices_ind]==0)[0]]
        return v, U_cp_v, U_nc_v
    elif flag=='out':
        U_cp_v = no_v_vertices_ind[ np.where(A[no_v_vertices_ind, v, ]>0)[0]]
        U_nc_v = no_v_vertices_ind[np.where(A[no_v_vertices_ind,v]==0)[0]]
        return v, U_cp_v, U_nc_v
    else:
        U_cp_v_in  = no_v_vertices_ind[ np.where(A[v, no_v_vertices_ind, ]>0)[0]]
        U_nc_v_in = no_v_vertices_ind[np.where(A[v, no_v_vertices_ind]==0)[0]]
        U_cp_v_out  = no_v_vertices_ind[ np.where(A[no_v_vertices_ind, v, ]>0)[0]]
        U_nc_v_out = no_v_vertices_ind[np.where(A[no_v_vertices_ind, v, ]==0)[0]]
        return v, U_cp_v_in, U_nc_v_in, U_cp_v_out, U_nc_v_out


# adaptive rewiring + random rewiring
def rewire_random(r,p,tau,A,v,U_cp_v,U_nc_v,flag,flag_alg):
    if r < p: # random
        j_minus = np.random.choice(U_cp_v)
        j_add = np.random.choice(U_nc_v)
    else: # diffusion
        if flag_alg == 'consensus':
            H = -compute_consensus_kernel(A, tau)
        elif flag_alg == 'advection':
            H = -compute_advection_kernel(A, tau)
        
        if flag == 'out':
            j_minus = U_cp_v[np.argmax(H[U_cp_v,v])]
            j_add = U_nc_v[np.argmin(H[U_nc_v,v])]
        else:
            j_minus = U_cp_v[np.argmax(H[v,U_cp_v])]
            j_add = U_nc_v[np.argmin(H[v,U_nc_v])]
    return j_add, j_minus



# adaptive rewiring + distance-based rewiring
def rewire_distance(r,p,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg):
    if r < p: # distance      
        j_minus = U_cp_v[np.argmax(D[U_cp_v,v])]
        j_add = U_nc_v[np.argmin(D[U_nc_v,v])]
    else: # diffusion
        if flag_alg == 'consensus':
            H = -compute_consensus_kernel(A, tau)
        elif flag_alg == 'advection':
            H = -compute_advection_kernel(A, tau)
        
        if flag == 'out':
            j_minus = U_cp_v[np.argmax(H[U_cp_v,v])]
            j_add = U_nc_v[np.argmin(H[U_nc_v,v])]
        else:
            j_minus = U_cp_v[np.argmax(H[v,U_cp_v])]
            j_add = U_nc_v[np.argmin(H[v,U_nc_v])]
    return j_add, j_minus






# adaptive rewiring + Boost-based rewiring
def rewire_boost(tau,A,v,U_cp_v,U_nc_v,flag,flag_alg, pBoost,boost_node, booster): 
    A_1 = A.copy()

    #Manipulate Adjacency Matrix to introduce a Boost (only in the node where calcium wave takes place given a probability)
    #booster = 1

    idx1 = np.nonzero(A_1[:,boost_node])
    A_1[idx1[0],boost_node] = A_1[idx1[0],boost_node]+booster
            

    # Weight renormalization because boost has taken place
    n_vertices = A_1.shape[0]
    edges = int(np.round(2 * np.log(n_vertices) * (n_vertices - 1), decimals=0))
    norm_factor = edges / np.sum(A_1)
    A_1 = A_1 * norm_factor
        
    
    # Advection vs Consensus
    if flag_alg == 'consensus':
        H = -compute_consensus_kernel(A_1, tau)
    elif flag_alg == 'advection':
        H = -compute_advection_kernel(A_1, tau)
    
    if flag == 'out':
        ########## Code fixing J's bug #################
        states_cp = H[U_cp_v,v]
        state_min = np.max(states_cp)
        states_nc = H[U_nc_v,v]
        state_max = np.min(states_nc)
        tie_cp = np.where(states_cp==state_min)[0]
        tie_nc = np.where(states_nc==state_max)[0]
        
        if len(tie_nc)>1:
            ind_add = np.random.choice(tie_nc)
        else:
            ind_add = tie_nc[0]
        
        if len(tie_cp)>1:
            ind_minus = np.random.choice(tie_cp)
        else:
            ind_minus = tie_cp[0]
            
        #tie_nc_flag = int(len(tie_nc)>1)
        #tie_cp_flag = int(len(tie_cp)>1)
        # adaptive rewiring
        j_minus = U_cp_v[ind_minus]
        j_add = U_nc_v[ind_add]
        #################################################
    
        # j_minus = U_cp_v[np.argmax(H[U_cp_v,v])]
        # j_add = U_nc_v[np.argmin(H[U_nc_v,v])]
    else:
        ########## Code fixing Jias bug #################
        states_cp = H[v,U_cp_v]
        state_min = np.max(states_cp)
        states_nc = H[v,U_nc_v]
        state_max = np.min(states_nc)
        tie_cp = np.where(states_cp==state_min)[0]
        tie_nc = np.where(states_nc==state_max)[0]
            
        if len(tie_nc)>1:
            ind_add = np.random.choice(tie_nc)
        else:
            ind_add = tie_nc[0]
        
        if len(tie_cp)>1:
            ind_minus = np.random.choice(tie_cp)
        else:
            ind_minus = tie_cp[0]
            
        #tie_nc_flag = int(len(tie_nc)>1)
        #tie_cp_flag = int(len(tie_cp)>1)
        # adaptive rewiring
        j_minus = U_cp_v[ind_minus]
        j_add = U_nc_v[ind_add]
        #################################################
        
        # j_minus = U_cp_v[np.argmax(H[v,U_cp_v])]
        # j_add = U_nc_v[np.argmin(H[v,U_nc_v])]
    return j_add, j_minus   #, H





# adaptive rewiring + boost-based rewiring + distance-based rewiring
def rewire_distance_boost(r,pDist,tau,D,A,v,U_cp_v,U_nc_v,flag,flag_alg, pBoost,boost_node, booster):
    
    if r < pDist: # distance    
        #print('Distance')
        j_minus = U_cp_v[np.argmax(D[U_cp_v,v])]
        j_add = U_nc_v[np.argmin(D[U_nc_v,v])]
    
    elif r < (pDist+pBoost): # Calcium wave
        #print('Boost')
        j_add, j_minus = rewire_boost(tau,A,v,U_cp_v,U_nc_v,flag,flag_alg, pBoost,boost_node, booster)

    else: # diffusion
        #print('Diffusion')
        if flag_alg == 'consensus':
            H = -compute_consensus_kernel(A, tau)
        elif flag_alg == 'advection':
            H = -compute_advection_kernel(A, tau)
        
        if flag == 'out':
            ########## Code fixing J's bug #################
            states_cp = H[U_cp_v,v]
            state_min = np.max(states_cp)
            states_nc = H[U_nc_v,v]
            state_max = np.min(states_nc)
            tie_cp = np.where(states_cp==state_min)[0]
            tie_nc = np.where(states_nc==state_max)[0]
            
            if len(tie_nc)>1:
                ind_add = np.random.choice(tie_nc)
            else:
                ind_add = tie_nc[0]
            
            if len(tie_cp)>1:
                ind_minus = np.random.choice(tie_cp)
            else:
                ind_minus = tie_cp[0]
            #tie_nc_flag = int(len(tie_nc)>1)
            #tie_cp_flag = int(len(tie_cp)>1)
            # adaptive rewiring
            j_minus = U_cp_v[ind_minus]
            j_add = U_nc_v[ind_add]
            #################################################
            
            # j_minus = U_cp_v[np.argmax(H[U_cp_v,v])]
            # j_add = U_nc_v[np.argmin(H[U_nc_v,v])]
            
        else:
            ########## Code fixing J's bug #################
            states_cp = H[v,U_cp_v]
            state_min = np.max(states_cp)
            states_nc = H[v,U_nc_v]
            state_max = np.min(states_nc)
            tie_cp = np.where(states_cp==state_min)[0]
            tie_nc = np.where(states_nc==state_max)[0]
            
            if len(tie_nc)>1:
                ind_add = np.random.choice(tie_nc)
            else:
                ind_add = tie_nc[0]
            
            if len(tie_cp)>1:
                ind_minus = np.random.choice(tie_cp)
            else:
                ind_minus = tie_cp[0]
            #tie_nc_flag = int(len(tie_nc)>1)
            #tie_cp_flag = int(len(tie_cp)>1)
            # adaptive rewiring
            j_minus = U_cp_v[ind_minus]
            j_add = U_nc_v[ind_add]
            #################################################
            
            # j_minus = U_cp_v[np.argmax(H[v,U_cp_v])]
            # j_add = U_nc_v[np.argmin(H[v,U_nc_v])]
        
    return j_add, j_minus   #, H





















# rewire out-degree based on rewiring principles
def rewire_out_principles(r,p,q,D,W,tau,A,v,U_cp_v,U_nc_v,flag_alg):
    if r < p: # distance      
        j_minus = U_cp_v[np.argmax(D[U_cp_v,v])]
        j_add = U_nc_v[np.argmin(D[U_nc_v,v])]
    elif r< p+q: # wave
        j_minus = U_cp_v[np.argmax(W[U_cp_v,v])]
        j_add = U_nc_v[np.argmin(W[U_nc_v,v])]
    else: # diffusion
        if flag_alg == 'consensus':
            H = -compute_consensus_kernel(A, tau)
        elif flag_alg == 'advection':
            H = -compute_advection_kernel(A, tau)
        states_cp = H[U_cp_v,v]
        state_min = np.max(states_cp)
        states_nc = H[U_nc_v,v]
        state_max = np.min(states_nc)
        tie_cp = np.where(states_cp==state_min)[0]
        tie_nc = np.where(states_nc==state_max)[0]
        
        if len(tie_nc)>1:
            ind_add = np.random.choice(tie_nc)
        else:
            ind_add = tie_nc[0]
        
        if len(tie_cp)>1:
            ind_minus = np.random.choice(tie_cp)
        else:
            ind_minus = tie_cp[0]
        #tie_nc_flag = int(len(tie_nc)>1)
        #tie_cp_flag = int(len(tie_cp)>1)
        # adaptive rewiring
        j_minus = U_cp_v[ind_minus]
        j_add = U_nc_v[ind_add]
    return j_add, j_minus


# rewire in-degree based on rewiring principles
def rewire_in_principles(r,p,q,D,W,tau,A,v,U_cp_v,U_nc_v,flag_alg):
    if r < p: # distance        
        j_minus = U_cp_v[np.argmax(D[v,U_cp_v])]
        j_add = U_nc_v[np.argmin(D[v,U_nc_v])]
    elif r< p+q: # wave
        j_minus = U_cp_v[np.argmax(W[v,U_cp_v])]
        j_add = U_nc_v[np.argmin(W[v,U_nc_v])]
    else: # diffusion
        if flag_alg == 'consensus':
            H = -compute_consensus_kernel(A, tau)
        elif flag_alg == 'advection':
            H = -compute_advection_kernel(A, tau)
        states_cp = H[v,U_cp_v]
        state_min = np.max(states_cp)
        states_nc = H[v,U_nc_v]
        state_max = np.min(states_nc)
        tie_cp = np.where(states_cp==state_min)[0]
        tie_nc = np.where(states_nc==state_max)[0]
        
        if len(tie_nc)>1:
            ind_add = np.random.choice(tie_nc)
        else:
            ind_add = tie_nc[0]
        
        if len(tie_cp)>1:
            ind_minus = np.random.choice(tie_cp)
        else:
            ind_minus = tie_cp[0]
        #tie_nc_flag = int(len(tie_nc)>1)
        #tie_cp_flag = int(len(tie_cp)>1)
        # adaptive rewiring
        j_minus = U_cp_v[ind_minus]
        j_add = U_nc_v[ind_add]
    return j_add, j_minus
