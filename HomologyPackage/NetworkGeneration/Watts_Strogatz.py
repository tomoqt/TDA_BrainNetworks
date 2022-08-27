# -*- coding: utf-8 -*-
import numpy as np

'''
def closest_nodes(n,dist,upper=False,non_diag=True):
  
    nodes=np.arange(n)
    col,row=np.meshgrid(nodes,nodes)

    closest=np.logical_or(np.abs(col-row)<=dist,np.abs(col-row)>=n-dist)
    if non_diag:
        closest=np.logical_and(closest, col!=row)
    if not upper:
        return  closest
  
    return np.logical_and(closest, col>row)

def old_nodes(old_graph, upper=False):
    n=len(old_graph)
    nodes=np.arange(n)
    col,row=np.meshgrid(nodes,nodes)
    old_nodes=old_graph>0
    if upper:
        return np.logical_and(col>row,old_nodes)
    return old_nodes

def new_nodes(old_graph,new_graph, upper=False):
    n=len(old_graph)
    nodes=np.arange(n)
    col,row=np.meshgrid(nodes,nodes)
    old=old_nodes(old_graph)
    new_nodes=np.logical_and(new_graph>0,np.logical_not(old))
    if upper:
        return np.logical_and(col>row,new_nodes)
    return new_nodes

def possible_nodes(new_graph, upper=False):
    n=len(new_graph)
    nodes=np.arange(n)
    col,row=np.meshgrid(nodes,nodes)
    possible_nodes= new_graph==0
    if upper:
        return np.logical_and(possible_nodes,col>row)
    return possible_nodes

def create_random_small_world_link(old_graph, new_graph,p,rng=0):
    # set seed
    if rng==0:
        rng=np.random.default_rng()

    n=len(old_graph)
    nodes=np.arange(n)
    cols,rows=np.meshgrid(nodes,nodes)

    olds=old_nodes(old_graph, upper=True)
    news=new_nodes(old_graph, new_graph, upper=True)
    possibles=possible_nodes(new_graph, upper=True)
    #shuffled is a np array where we save all the nodes of the newgraph after a shuffle
    shuffled=np.zeros(np.shape(news))
    shuffled[olds]=1
    shuffled[news]=rng.uniform(0,1,len(shuffled[news]))
    shuffled[shuffled>p]=1

    #now i remain with all the nodes I need to change
    for row in range(n):
        for col in range(n):
            if shuffled[row,col]>0 and shuffled[row,col]<1:

                if np.any(possibles[row]):

                    appo_index=rng.integers(0,len(shuffled[possibles[row]]))
                    index=cols[row,possibles[row]][appo_index]

                    shuffled[row,index]=1
                    shuffled[row,col]=0
                    possibles=possible_nodes(shuffled, upper=True)


    return shuffled+np.transpose(shuffled)



from ring_lattice import ring_lattice

def Watts_Strogatz(n,p, seed=np.nan):
    max_edge=n*(n-1)/2

    # preallocate output
    corr_matrix=np.zeros([n,n])

    # generate random correlation

    # set seed
    if  np.isnan(seed):
        rng=np.random.default_rng()
        seed=rng.integers(1000)        
    rng=np.random.default_rng(seed)

    a=1
    G1=np.zeros([n,n])
    G_adjacency=np.zeros([n,n])

    for k in range(int(np.floor(n/2))):
        #generate closest
        G1=G_adjacency
        G2=closest_nodes(n,dist=k+1).astype(int)#,upper=True)

    #shuffle
        G_adjacency=create_random_small_world_link(G1,G2,p,rng)

        G_appo=np.triu(G_adjacency)-np.triu(G1)

        length=sum(sum(G_appo==1))
        G_appo[G_appo==1]=rng.uniform(a-length/max_edge,a,length)
        a=a-length/max_edge

        corr_matrix+=G_appo

    return corr_matrix+np.transpose(corr_matrix)
'''

from .ring_lattice import ring_lattice
def Watts_Strogatz(n,p,seed=np.nan):
    # preallocate output
    

    # generate random correlation

    # set seed
    if  np.isnan(seed):
        rng=np.random.default_rng()
        seed=rng.integers(1000)        
    rng=np.random.default_rng(seed)
    
    # preallocate output
    corr_matrix=ring_lattice(n,seed=seed)
    corr_matrix=np.triu(corr_matrix)
    
    #find_rewiring nodes
    rewiring_nodes=rng.uniform(0,1,[n,int(np.ceil(n/2))-2])<p
    for dist in np.arange(np.shape(rewiring_nodes)[1])+1:
        for node in np.arange(n):
            if node+dist>=n:
                d=dist-n
            else:
                d=dist
                
            if rewiring_nodes[node][dist-1]:
                arrival_node_dist=rng.integers(dist+1,n-dist-1)
                arrival_node=node+arrival_node_dist
                
                if arrival_node>=n:
                    arrival_node=arrival_node-n
                    if d==dist:
                        corr_matrix[node][node+d],corr_matrix[arrival_node][node]=corr_matrix[arrival_node][node],corr_matrix[node][node+d]
                    else:
                        corr_matrix[node+d][node],corr_matrix[arrival_node][node]=corr_matrix[arrival_node][node],corr_matrix[node+d][node]
                
                else:
                    if d==dist:
                        corr_matrix[node][node+d],corr_matrix[node][arrival_node]=corr_matrix[node][arrival_node],corr_matrix[node][node+d]
                    else:
                        corr_matrix[node+d][node],corr_matrix[node][arrival_node]=corr_matrix[node][arrival_node],corr_matrix[node+d][node]
                        
                    
                
    return corr_matrix+corr_matrix.transpose()
                

#Watts_Strogatz(10,0.2)