'''
From correlation martix to graph:

This modules has all the function necessaries to get a graph given a correlation matrix and a threshold eps

The first funcion gives us a not connected graph, with the second one we are assured to only one connected component 
'''


import networkx as nx
#import math
from networkx import tree
import numpy as np
import scipy.special as ss

def normalize(corr_mat):
    nnodes=int(len(corr_mat))
    np.fill_diagonal(corr_mat,0)
    corr_mat=np.triu(corr_mat)
    corr_mat_flatten=corr_mat.reshape([nnodes**2])
    nlinks=int(nnodes*(nnodes-1)/2)
    indices=np.argsort(corr_mat_flatten)
    values=np.linspace(1/(nlinks),1,nlinks)
    values=np.concatenate((np.zeros(int(nnodes**2-nlinks)),values))
    corr_mat_flatten[indices]=values
    
    return corr_mat_flatten.reshape([nnodes,nnodes])+   corr_mat_flatten.reshape([nnodes,nnodes]).transpose()  
    
    


def adj_matrix_notconnected(corr_mat, eps):#this function it is not used, but maybe it will be useful in the future

    normalized_corr_mat=normalize(corr_mat)
    adj_matrix=normalized_corr_mat>1-eps

    return adj_matrix.astype(int) 



def adj_matrix_connected(corr_matrix,sparsity_value):
    #non funziona attualmente, non so perchè
    """
    given the correlation matrix and the expected sparsity coefficient it can 
    happen that the corresponding thresholded matrix results in a disconnected graph
    here we force the graph to be fully connected by the computation of the minimum
    spanning tree and adding the required edges in order to have a unique connected component 
    """
    if sparsity_value == 1.0:
        adj_matrix=np.ones(corr_matrix.shape)
        np.fill_diagonal(adj_matrix,0)
        return adj_matrix
        
    #aggiungere l'errore se lo sparsity value è troppo basso
    
    corr_matrix =abs(corr_matrix)

    max_num_edges = ss.comb(corr_matrix.shape[0],2)
    num_edges = int(max_num_edges*sparsity_value)
    
    num_regions=corr_matrix.shape[0]
    #total number of regions in the graph


    totalgraph=nx.from_numpy_matrix(1-abs(corr_matrix))
    #extraction of a complete graph having has weight 1-abs(correlation)
    #we need to take 1-abs since the mst is taking the minimum weight graph and we want the most correlated edges to be there
    
    MST=nx.to_numpy_array(tree.minimum_spanning_tree(totalgraph).to_undirected())
    MST_adj_mat=MST #inutile
    MST_adj_mat[MST>0]==1 #inutile
    MST_adj_mat=np.triu(MST_adj_mat) #put zeros in the inferior triangular matrix
    
    #put zeros in the diagonal of the corr matrix
    for i in range(num_regions):
        corr_matrix[i,i]=0
    
    values_corr=abs(np.triu(corr_matrix)) #forse abs inutile
    
    cor_wo_MST=values_corr[np.triu(MST_adj_mat)==0]
    #we do not consider the correlation values which do not involve edges that are already in the MST
        #we consider the correlation values which do not involve edges that are already in the MST ???
    values=list(cor_wo_MST.flatten())
    values.sort(reverse=True)
    
    #we select the maximum value of correlation to have the expected num of edges - num of edges in the mst (num regions-1)
    # values must have at least 1 element: AT LEAST 1 LINK MORE THAN  THE  SPANNING TREE
    value_thresh=values[min(num_edges-(num_regions-1),len(values))-1] #-1 index start at 0
    #notice, we can't choose a sparsity value resulting in a nuber of edges lower than then number of edges in the spanning tree (sparsity value tree  = 2*(n-1)/(n-1)*n = 1*2/n)
    adj_matrix=np.zeros(corr_matrix.shape) 
    
    #we put an edge if the value of correlation is higher than the found threshold or if the edges is required by the mst
    adj_matrix[values_corr>=value_thresh]=1
    adj_matrix[MST_adj_mat!=0]=1
    
    adj_matrix=(adj_matrix)+np.transpose((adj_matrix)) #simmetry of the adj matrix
    
    return adj_matrix


def adjacency_matrix(corr_mat,eps,method='threshold'):
    #verify corr_mat is a sqare array
    
    cond= len(np.shape(corr_mat))==2 and np.shape(corr_mat)[0]==np.shape(corr_mat)[1]
    if cond:
        if method=='sparsity_connected':
            if np.floor(eps*len(corr_mat))>=len(corr_mat)-1:
                return adj_matrix_connected(corr_mat, eps)
            else:
                print('sparsity value is too low')
        elif method=='sparsity_notconnected':
            return adj_matrix_connected(corr_mat, eps)
        elif method=='threshold':
            np.fill_diagonal(corr_mat, 0)
            mask=corr_mat>eps
            adj_mat=np.zeros(np.shape(corr_mat))
            adj_mat[mask]=1
            return adj_mat
        else:
            method=input('try methods: \n sparsity_connected \n sparsity_notconnected \n threshold')
            return adjacency_matrix(corr_mat,eps,method=method)
    else:
        print('corr_mat is not a 2-dimensional square array')
              
        

