import numpy as np
import pickle

#import sys
#sys.path.append(".")
#from scripts import rewire, indices, directed_modularity, initialization_params
import indices, directed_modularity

import os
os.chdir('../scripts')

f = open('./Output/boost_nor.pckl', 'rb')
A_nor = pickle.load(f)
f.close()

f = open('./Output/boosted_nodes.pckl', 'rb')
boosted_nodes = pickle.load(f)
f.close()

os.chdir('../Data_analysis')

# METRICS

#Connectivity and path lengths
path_length = {}; connected_pairs = {}

for tau in A_nor.keys():
    path_length[tau] = {}
    connected_pairs[tau] = {}
    for p_in in A_nor[tau].keys():
        path_length[tau][p_in] = {}
        connected_pairs[tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            path_length[tau][p_in][pDist] = {}
            connected_pairs[tau][p_in][pDist] = {} 
            for pBoost in A_nor[tau][p_in][pDist].keys():
                path_length[tau][p_in][pDist][pBoost] = {}
                connected_pairs[tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    path_length[tau][p_in][pDist][pBoost][k] = {}
                    connected_pairs[tau][p_in][pDist][pBoost][k] = {}
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                        connected_pairs[tau][p_in][pDist][pBoost][k][r], path_length[tau][p_in][pDist][pBoost][k][r] = indices.path_length_connectedness(A_nor[tau][p_in][pDist][pBoost][k][r])


fname = './Metrics/all_cp.pckl'
f = open(fname, 'wb')
pickle.dump(connected_pairs, f)
f.close()

fname = './Metrics/all_pl.pckl'
f = open(fname, 'wb')
pickle.dump(path_length, f)
f.close()



#Number of convergent and divergent hubs

binary_flag = True; thresh = 15 

Hubs = {}
Hubs['in'] = {}; Hubs['out'] = {}

for tau in A_nor.keys():
    Hubs['in'][tau] = {}
    Hubs['out'][tau] = {}
    for p_in in A_nor[tau].keys():
        Hubs['in'][tau][p_in] = {}
        Hubs['out'][tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            Hubs['in'][tau][p_in][pDist] = {}
            Hubs['out'][tau][p_in][pDist] = {} 
            for pBoost in A_nor[tau][p_in][pDist].keys():
                Hubs['in'][tau][p_in][pDist][pBoost] = {}
                Hubs['out'][tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    Hubs['in'][tau][p_in][pDist][pBoost][k] = {}
                    Hubs['out'][tau][p_in][pDist][pBoost][k] = {}
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                        Hubs['in'][tau][p_in][pDist][pBoost][k][r] = indices.hub_number(A_nor[tau][p_in][pDist][pBoost][k][r],thresh,binary_flag, axisUsed=1)
                        Hubs['out'][tau][p_in][pDist][pBoost][k][r] = indices.hub_number(A_nor[tau][p_in][pDist][pBoost][k][r],thresh,binary_flag, axisUsed=0)


fname = './Metrics/all_hubs_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(Hubs, f)
f.close()




#Obtain node numbers of convergent and divergent hubs

Hubs_nodeNo = {}
Hubs_nodeNo['in'] = {}; Hubs_nodeNo['out'] = {}

for tau in A_nor.keys():
    Hubs_nodeNo['in'][tau] = {}
    Hubs_nodeNo['out'][tau] = {}
    for p_in in A_nor[tau].keys():
        Hubs_nodeNo['in'][tau][p_in] = {}
        Hubs_nodeNo['out'][tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            Hubs_nodeNo['in'][tau][p_in][pDist] = {}
            Hubs_nodeNo['out'][tau][p_in][pDist] = {} 
            for pBoost in A_nor[tau][p_in][pDist].keys():
                Hubs_nodeNo['in'][tau][p_in][pDist][pBoost] = {}
                Hubs_nodeNo['out'][tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    Hubs_nodeNo['in'][tau][p_in][pDist][pBoost][k] = {}
                    Hubs_nodeNo['out'][tau][p_in][pDist][pBoost][k] = {}
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                                
                        Topo_A = 1.0*(A_nor[tau][p_in][pDist][pBoost][k][r]>0)
                        deg_in = np.sum(Topo_A,1)
                        deg_out = np.sum(Topo_A,0)    
                        Con = np.where((deg_in>=thresh)&(deg_out>0))[0]
                        Div = np.where((deg_out>=thresh)&(deg_in>0))[0]
                                
                        Hubs_nodeNo['in'][tau][p_in][pDist][pBoost][k][r] = Con
                        Hubs_nodeNo['out'][tau][p_in][pDist][pBoost][k][r] = Div


fname = './Metrics/all_hubs_nodeNo_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(Hubs_nodeNo, f)
f.close()




#Obtain each node's in-degree and out-degree

Degree = {}
Degree['in'] = {}; Degree['out'] = {}

for tau in A_nor.keys():
    Degree['in'][tau] = {}
    Degree['out'][tau] = {}
    for p_in in A_nor[tau].keys():
        Degree['in'][tau][p_in] = {}
        Degree['out'][tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            Degree['in'][tau][p_in][pDist] = {}
            Degree['out'][tau][p_in][pDist] = {} 
            for pBoost in A_nor[tau][p_in][pDist].keys():
                Degree['in'][tau][p_in][pDist][pBoost] = {}
                Degree['out'][tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    Degree['in'][tau][p_in][pDist][pBoost][k] = {}
                    Degree['out'][tau][p_in][pDist][pBoost][k] = {}
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                                                                    
                        Topo_A = 1.0*(A_nor[tau][p_in][pDist][pBoost][k][r]>0)
                        deg_in = np.sum(Topo_A,1)
                        deg_out = np.sum(Topo_A,0)
                        nodes = np.arange(0,100)
                             
                        Degree['in'][tau][p_in][pDist][pBoost][k][r] = np.stack((nodes,deg_in), axis=0)
                        Degree['out'][tau][p_in][pDist][pBoost][k][r] = np.stack((nodes,deg_out), axis=0)
                                


fname = './Metrics/InOut_degree.pckl'
f = open(fname, 'wb')
pickle.dump(Degree, f)
f.close()




#Get cd-units

f = open('./Metrics/all_cp.pckl', 'rb')
connected_pairs = pickle.load(f)
f.close()

cd_pairs = {}

for tau in A_nor.keys():
    cd_pairs[tau] = {}
    for p_in in A_nor[tau].keys():
        cd_pairs[tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            cd_pairs[tau][p_in][pDist] = {}
            for pBoost in A_nor[tau][p_in][pDist].keys():
                cd_pairs[tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    cd_pairs[tau][p_in][pDist][pBoost][k] = {}
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                        cd_pairs[tau][p_in][pDist][pBoost][k][r] = indices.cd_pairs(A_nor[tau][p_in][pDist][pBoost][k][r],connected_pairs[tau][p_in][pDist][pBoost][k][r], thresh)


fname = './Metrics/all_cd_pairs_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(cd_pairs, f)
f.close()




#Size and density of the intermediate subgraphs

f = open('./Metrics/all_cd_pairs_'+str(thresh)+'.pckl', 'rb')
connected_cd_pairs = pickle.load(f)
f.close()


interm_size = {}; interm_density = {}; periph_density = {}


p_in = 0.5
for tau in A_nor.keys():
    interm_size[tau] = {}
    interm_density[tau] = {}
    periph_density[tau] = {}
    for p_in in A_nor[tau].keys():
        interm_size[tau][p_in] = {}
        interm_density[tau][p_in] = {}
        periph_density[tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            interm_size[tau][p_in][pDist] = {}
            interm_density[tau][p_in][pDist] = {} 
            periph_density[tau][p_in][pDist] = {} 
            for pBoost in A_nor[tau][p_in][pDist].keys():
                interm_size[tau][p_in][pDist][pBoost] = {}
                interm_density[tau][p_in][pDist][pBoost] = {}
                periph_density[tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    interm_size[tau][p_in][pDist][pBoost][k] = {}
                    interm_density[tau][p_in][pDist][pBoost][k] = {}
                    periph_density[tau][p_in][pDist][pBoost][k] = {}
                            
                    for ind,r in enumerate(A_nor[tau][p_in][pDist][pBoost][k].keys()):
                        if connected_pairs[tau][p_in][pDist][pBoost][k][r].shape[1]>0:
                            interm_size[tau][p_in][pDist][pBoost][k][r],interm_density[tau][p_in][pDist][pBoost][k][r],periph_density[tau][p_in][pDist][pBoost][k][r] = indices.intermediate_subgraphs(A_nor[tau][p_in][pDist][pBoost][k][r], connected_cd_pairs[tau][p_in][pDist][pBoost][k][r])
                        else:
                            interm_size[tau][p_in][pDist][pBoost][k][r] = np.nan
                            interm_density[tau][p_in][pDist][pBoost][k][r] = np.nan
                            periph_density[tau][p_in][pDist][pBoost][k][r] = np.nan
                            


fname = './Metrics/all_interm_size_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(interm_size, f)
f.close()

fname = './Metrics/all_interm_density_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(interm_density, f)
f.close()

fname = './Metrics/all_periph_density_'+str(thresh)+'.pckl'
f = open(fname, 'wb')
pickle.dump(periph_density, f)
f.close()



#Modularity

modularity = {}
mean_modularity_k= {}
mean_modularity_total= {}

p_in = 0.5


import statistics

for tau in A_nor.keys():
    modularity[tau] = {}
    mean_modularity_k[tau] = {}
    mean_modularity_total[tau] = {}
    for p_in in A_nor[tau].keys():
        modularity[tau][p_in] = {}
        mean_modularity_k[tau][p_in] = {}
        mean_modularity_total[tau][p_in] = {}
        for pDist in A_nor[tau][p_in].keys():
            modularity[tau][p_in][pDist] = {}
            mean_modularity_k[tau][p_in][pDist] = {}
            mean_modularity_total[tau][p_in][pDist] = {}
            for pBoost in A_nor[tau][p_in][pDist].keys():
                modularity[tau][p_in][pDist][pBoost] = {}
                mean_modularity_k[tau][p_in][pDist][pBoost] = {}
                mean_modularity_total[tau][p_in][pDist][pBoost] = {}
                for k in A_nor[tau][p_in][pDist][pBoost].keys():
                    modularity[tau][p_in][pDist][pBoost][k] = {}
                            
                            
                    sum_vec = 0
                    for r in A_nor[tau][p_in][pDist][pBoost][k].keys():
                        (modularity[tau][p_in][pDist][pBoost][k][r], communitiesDict) = directed_modularity.getModularityIndex(A_nor[tau][p_in][pDist][pBoost][k][r])

                        sum_vec = sum_vec+np.matrix(modularity[tau][p_in][pDist][pBoost][k][r])
                                    
                                    
                    mean_modularity_k[tau][p_in][pDist][pBoost][k] = sum_vec/len(A_nor[tau][p_in][pDist][pBoost][k])
                            

                mean_modularity_total[tau][p_in][pDist][pBoost] = np.mean(list(mean_modularity_k[tau][p_in][pDist][pBoost].values()))



fname = './Metrics/mean_modularity.pckl'
f = open(fname, 'wb')
pickle.dump(mean_modularity_total, f)
f.close()

fname = './Metrics/all_modularity.pckl'
f = open(fname, 'wb')
pickle.dump([mean_modularity_total,mean_modularity_k,modularity], f)
f.close()