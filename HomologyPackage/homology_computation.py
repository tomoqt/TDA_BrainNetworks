'''
HOMOLOGY_COMPUTATION

In this module we implemented all the function used to compute topological
properties of networks

First we implemented our custom vietoris-rips persistance diagrams, using GIOTTO TDA.
All the other function are based on this first function


'''

import numpy as np
import networkx as nx
from .filtration import *
from gtda.homology import VietorisRipsPersistence
import cechmate as cm
from timeit import default_timer as timer




def custom_vietoris_persistance (corr_mats,costs,homology_dimensions,plot=False):
    '''
    input:
    corr_mats:  correlation matrices of the graphs we want to compute the diagram
    costs:      list of the costs where the filtration is computed
    homology_dimensions: list of the dimensions of the homology to compute
    plot:       boolean value, if true this function plot the persistence graph 
  
    output:       
    diagrams: list of persistance diagrams. Each element of the list 
              is a list of points in the plane with associated dimansion
    '''  

    if isinstance(homology_dimensions, int):
        homology_dimensions=list(range(1,homology_dimensions+1))

    if not isinstance(corr_mats,list):
        corr_mats=[corr_mats]

    #generate filtrations
    graph_filtrations = [generate_graph_sequence(g,costs) for g in corr_mats]
    #generate adjacency matrices
    adjacency_matrices = [from_filtration_to_adjacency(graph_filtration) for graph_filtration in graph_filtrations]
    # Instantiate topological transformer
    VR = VietorisRipsPersistence(metric="precomputed",homology_dimensions=homology_dimensions)

    # Compute persistence diagrams corresponding to each graph in X
    diagrams = VR.fit_transform(adjacency_matrices)

    if (plot):
        VR.plot(diagrams, sample=0)
  
    return diagrams,VR



def custom_filtration_persistance(corr_mats,costs,homology_dimensions):
    '''
    input:
    corr_mats:  correlation matrices of the graphs we want to compute the diagram
    costs:      list of the costs of which we extract a graph
    homology_dimensions: list of the dimensions of the homology to compute
    plot:       boolean value, if true this function plot the persistence graph 
  
    output:       
    diagrams: list of list of persistance diagrams. 
    Each element of the outer list represent a brain network,
    each element of the inner list 
    is a list of points in the plane with associated dimansion
    '''  

    if isinstance(homology_dimensions, int):
        homology_dimensions=list(range(1,homology_dimensions+1))

    if not isinstance(corr_mats,list):
        corr_mats=[corr_mats]

    # graph list is a list of all the graph for each brain
    graph_list = [generate_graph_sequence(g,costs) for g in corr_mats]
    #print(graph_list)
    print("graph list fatta")
    #filtrations
    #filtrations = [[custom_filtration(G,homology_dimensions) for G in brain] for brain in graph_list]
    filtrations=[]
    for i,brain in enumerate(graph_list):
        print("brain" , i/len(graph_list)*100)
        filtrations.append([])
        for j,G in enumerate(brain):
            print(j/len(brain)*100)
            start=timer()
            filtrations[i].append(custom_filtration(G,homology_dimensions))
            end=timer()
            print(start-end)
        
    #print(len(filtrations[0][1]))
    #print(len(filtrations[0][0]))

    dgms = [[cm.phat_diagrams(filtration, show_inf = True) for filtration in brain ] for brain in filtrations]

    return (dgms,filtrations)


'''
BETTI's NUMBERS
'''

def separate(diagrams):
    '''
    input:  
    dmg: persistence diagram from giotto 
    output: 
    H: list of different dimensions points in persistent diagrams
    '''
    H=[]

    for i,diagram in enumerate(diagrams):
    
        H.append([])
        #set of the dimensions of the points in the diagram
        dimensions=set(diagram[:,2])
    
        #for all the dimension i generate the list of the point of that exact dimension
        for j,dimension in enumerate(dimensions):
            H[i].append([])
            H[i][j]=[point for point in diagram if point[2]==dimension]
            H[i][j]=np.asarray(H[i][j])

            #if there is only one diagram there is no need of creating a list of list
            #if len(diagrams)==1:
            #  return H[0]

    return H


def betti_numbers (diagrams,n_filtration_steps): 
    '''
    input: 
    diagrams: persistence diagram from giotto 
    n_filtration_steps: it is the number of filtration steps
    output:
    betti_numbers: it is a list of matrix where each column represent a dimension 
        and the rows are the filtration steps
    '''
    separated_diagram=separate(diagrams)
    betti_numbers=[]
    for k,H in enumerate(separated_diagram):
        n_dimensions=len(H)   # number of different dimensions where we computed the homology
        betti_numbers.append(np.zeros(shape=(n_filtration_steps, n_dimensions)))

        for i in range(n_filtration_steps):
            for j in range(n_dimensions):
                born = H[j][:,0]<=i+1 #generator of the omology must be already born
                not_death = i+1 < H[j][:,1] #generator og the homology should not be already dead
                betti_numbers[k][i,j]=sum(np.logical_and(born, not_death))
  
    return betti_numbers


'''
EULER CHARACTERISTIC, evaulated by brute force
'''


def euler_char_cliques(adj):
    #qua stiamo facendo forza bruta, magari possiamo calcolarci imparando dalle cricche della filtrazione precedente
    cliques=list(nx.enumerate_all_cliques(nx.from_numpy_matrix(adj))) # this should give the same value as its analogue w/ the betti numbers
    x=0
    dim_cliques=np.asarray([len(clique)-1 for clique in cliques])
    x=sum((-1)**dim_cliques)
    # for i in range(len(cliques)):
    #  x += ((-1)**(len(cliques[i])-1)) 
    return x

def euler_entropy_cliques(adj):
    return np.log(np.abs(euler_char_cliques(adj)))


'''
(Partial) EULER CHARACTERISTIC, evaulated from betti's number 
'''

def euler_char_betti(betti):
    # problema, noi non ci siamo calcolati tutti tutti i numeri di betti, solo nelle dimensioni che abbiamo scelto, che se fa?
    x=0
    a = np.empty((len(betti)))
    a[::2] = 1
    a[1::2] = -1
    # a è un vettore con componenti (-1)^n
    x = np.sum(np.multiply(betti,a))
    return x

def euler_entropy_betti(betti):
    return np.log(np.abs(euler_char_betti(betti)))
















