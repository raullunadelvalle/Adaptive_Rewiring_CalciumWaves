# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:18:33 2022

@author: Raúl Luna
"""

import numpy as np
import random





def boost_node_selector(init_mat):
    
    number_nodes = init_mat.shape[0]
    
    deg_out = np.sum(init_mat > 0, axis=0, keepdims=False,)
    nodes_sending = np.where((deg_out > 0) & (deg_out < number_nodes - 1))
    
    # the node selected for boosts must have at least 1 out-link
    boost_node = -9999
    while boost_node not in nodes_sending[0]:  
        rand_idx = random.randrange(len(nodes_sending[0]))
        boost_node = [nodes_sending[0][rand_idx]]  
            
    #boost_nodex = [np.random.randint(0, number_nodes-1)]

    return boost_node








