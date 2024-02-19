# CREATE INITIAL RANDOM NETWORKS
import numpy as np
import pickle
import rewire

# network settings
n_vertices = 100 
edges = int(np.round(2 * np.log(n_vertices) * (n_vertices - 1), decimals=0))

# number of different configurations
K = 150

# Weight distribution: normal
vertices_coord = {}
initial_matrix = {}
for k in range(K):
    vertices_coord[k], initial_matrix[k] = rewire.initial_directed_network(n_vertices, edges, 
                                                                           weight_distribution= 'normal',
                                                                           mu=1.0,sig=0.25)

f = open('./Output/initials_nor.pckl', 'wb') 
pickle.dump([vertices_coord,initial_matrix], f)
f.close()