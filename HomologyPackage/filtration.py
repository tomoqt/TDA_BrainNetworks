### FILTRATION.PY
'''
In this module we define all the function about the filtrations.
First we work on the filtration induced by threshold eps, then on nighborhods
'''



import numpy as np
#from .from_corr_to_graph import adj_matrix_connected
#from .from_corr_to_graph import adj_matrix_notconncted
#import networkx as nx
#from timeit import default_timer as timer




###
# General use function that gives us a matrix in the format needed to compute the diagrams by giotto tda
###

def from_filtration_to_adjacency(graph_filtration, sparsity_values):
  #forse da modificare mettendo i costi invece di interi
    
  sparsity_values=np.sort(sparsity_values)
  # pre allocate output
  # we create an array with non feasible weigth
  adjacency_matrix=np.ones(np.shape(graph_filtration[0]))*(-1)

  # for each graph in the filtration if a new link appears, the weigth is the 
  # number of the filtration.
  for i,graph in enumerate(graph_filtration):
    cond=np.logical_and(adjacency_matrix==-1,graph==1)
    adjacency_matrix[cond]=sparsity_values[i]
  
  # diagonal entries need to be 0
  adjacency_matrix+= np.identity(np.shape(adjacency_matrix)[0])
  # all non feasible weight becomes infinite 
  adjacency_matrix[adjacency_matrix==-1]=np.inf
  return adjacency_matrix




###
# FILTRATION FROM THE THRESHOLDS(non mi sembra così utile)
###
'''
def generate_graph_sequence(corr_matrix, costs):

  # input:
  #   correlation matrix: is the correlation matrix where we need to extract 
  #                       using the function provided to us
  #   costs: is the list of cost where we decided to make the filtration
  # output:
  #   graph_sequence: it is a list of numpy matrices, each being an adjacency 
  #                   matrix of a graph of hte filtration
  
  graph_sequence=list()
  for eps in costs:
    #graph_sequence.append(adj_matrix_connected(corr_matrix,eps))
    graph_sequence.append(adj_matrix_notconncted(corr_matrix, eps))
  return graph_sequence

'''
###
# NEIGHBORHOOD (not done)
###


###
# New filtration
###

import igraph as ig
# forse sarà inutile

def custom_filtration(adjacency_matrix, homology_dimensions):
    max_dim=max(homology_dimensions)
    G=ig.Graph.Adjacency(adjacency_matrix)
    cliques=G.cliques(max=max_dim+2)
    cliques = sorted(cliques, key=len)
    filtration=[(list(clique),len(clique)-1) for clique in cliques]
    return filtration
  















