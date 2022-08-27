# Ring Lattice generator
import numpy as np

def ring_lattice(n, seed=np.nan):
    '''
    Input:
    n -> number of nodes
    seed -> seed for pseudo-random number generation

    Output:
    adj_matrix -> weighted adjacency matrix for ring lattice graph
    '''    
    #max number of indirect edges
    max_edge = int(n*(n-1)/2)
    # preallocate output
    adj_matrix=np.zeros([n,n])
    # set seed
    if  np.isnan(seed):
        rng=np.random.default_rng()
        seed=rng.integers(1000)        
    rng=np.random.default_rng(seed)
    nodes=np.arange(n)
    col,row=np.meshgrid(nodes,nodes)   
    triu=row<col
    max_dist=int(np.floor(n/2))
    a=1
    for k in range(1,max_dist+1):
        kdist=np.logical_or(np.abs(col-row)==k,np.abs(col-row)==n-k)
        kdist_triu=np.logical_and(kdist,triu)
        length=np.sum(kdist_triu)
        adj_matrix[kdist_triu]=rng.uniform(a-length/max_edge,a,length)
        a=a-length/max_edge
    return adj_matrix+np.transpose(adj_matrix) 
  