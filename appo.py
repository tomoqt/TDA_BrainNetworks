from HomologyPackage import *
from HomologyPackage import adj_matrix_notconnected
G1=NetworkGeneration.Erdos_Renyi(90)
NetworkGeneration.Watts_Strogatz(10,0.2)
adj_matrix_notconnected(G1,0.2)

G1=NetworkGeneration.Watts_Strogatz(20,0.2)
G=ig.Graph.Adjacency(adj_matrix_notconnected(G1,0.2))
ig.plot(G)
edges=[e.tuple for e in G.es]
simplex=SimplexTree()
for edge in edges:
    simplex.insert(list(edge))

simplex.dimension()
simplex.expansion(90)
simplex.dimension()
simplex.compute_persistence()
simplex.betti_numbers()

filtration=custom_filtration(adj_matrix_notconnected(G1,0.2),[1])
from gudhi import SimplexTree

simplex=SimplexTree()
for i in range(len(filtration)):
    simplex.insert(filtration[i][0],filtration[i][1])
    
    

from gtda.homology import VietorisRipsPersistence
VR = VietorisRipsPersistence(metric="precomputed",homology_dimensions=[1,2])
diagrams = VR.fit_transform([G1])
VR.plot(diagrams, sample=0)
#from gtda.homology import VietorisRipsPersistence
'''
import numpy as np
import networkx as nx
import math
from matplotlib import pyplot as plt
import pandas as pd
from homology_modules import *
import pandas as pd
import cechmate as cm

c_m=generate_small_world_correlation_matrix(1000,0.01) 
df=pd.DataFrame(cm)
df.to_csv('matricepermatlab.csv')
print(df)


#nx.draw(nx.from_numpy_matrix(adj_matrix_notconncted(G,0.1)))
#nx.draw(nx.from_numpy_matrix(adj_matrix_notconncted(G,0.3)))
def neigborhood(G,node):
    #node=0
    closest_neighbors=list(nx.all_neighbors(G,node))
    print(closest_neighbors)
    #closest_neighbors.append(node)
    print(closest_neighbors)
    nodes=closest_neighbors
    #nodes=set()
    #for n in closest_neighbors:
    #    nodes.update(nx.neighbors(G,n))
    #print(nodes)
    #nodes.remove(node)
    return G.subgraph(nodes)

#nx.draw(nx.from_numpy_matrix(adj_matrix_notconncted(G,0.5)))
G=nx.from_numpy_matrix(adj_matrix_notconncted(c_m,0.1))

nx.draw(neigborhood(G,68))


cm_brain = np.loadtxt('regional-differentiation-based-on-graph-nodal-statistics-for-functional-brain-connectivity-networks-characterization/DATA/cor_mat_coma/cor_mat_90_wave_scale_3_ts_Control_grey_matter_ROI_90_ts_Control_0.txt.txt',usecols=range(1,91), skiprows=1, dtype='float')
G_brain=nx.from_numpy_matrix(adj_matrix_notconncted(cm_brain,0.1))
nx.draw(neigborhood(G_brain,18))

nx.draw(G)
nx.draw(G_brain)
#ingenerale sembranoesserci buchi al centro instrogatz e in periferia a brain functional network
# ci si aspetta che per valori bassi non si riesca a distinguere il centro dalla periferia



#fisso un costo e vedo l'evoluzione dell'omologia al variare di beta
def count_inf(dgm,dim):
        if (len(dgm)<=dim):
                return 0
        out=dgm[dim]
        return len(out[out[:,1]==np.inf])

def homology(corr_mat,costs,dim):
        graph_list = generate_graph_sequence(corr_mat,costs) #lista di grafi per il costo per un singolo brain
        #nbhds_list=[[ neighborhood(g,n) for n in range(len(corr_mat))] for g in graph_list] #per ogni costo ho la lista dei sottografi dei vicini, ogni elemento corrisponde a un nodo
        #print(nbhds_list)
        filtrations=[]
        homology_dimensions=list(range(dim+1))
        for k,graph in enumerate(graph_list):
                print(f"calcolo neighborhoods corripondendenti al costo {k+1} di {len(graph_list)}")
                #filtrations.append()
                start=timer()
                #for i,nbhd in enumerate(nbhds):
                filtrations.append(custom_filtration(graph,homology_dimensions))
                end=timer()
                print(end-start)

        #print(filtrations)
        dgms = [cm.phat_diagrams(filtration, show_inf = True) for filtration in filtrations]
        infs=np.array([count_inf(dgm,dim) for dgm in dgms])
        return (dgms,infs,filtrations)



def betti(corr_mats, costs, dim):
        infs_list=[]
        for i,mat in enumerate(corr_mats):

                print(f"SIAMO ALLA MATRICE n.{i+1} di {len(corr_mats)}")
                (dgms,infs,filtrations)=homology(mat,costs,dim)
                infs_list.append(infs)
        return infs_list
                #with open(save_file, "ab") as f:
                #        pickle.dump((dgms,infs,filtrations), f, protocol=pickle.HIGHEST_PROTOCOL)

#p=0.02 maximize small world propemsity


small_worlds1=[generate_small_world_correlation_matrix(90,0.01*i, i+2) for i in range(100)]
betti_ns=betti(small_worlds1, [0.1], 1)

# questo è un grafico da fare
small_worlds1=[[generate_small_world_correlation_matrix(90,0.01*i, k) for i in range(100)] for k in range(100)]
betti_ns1=[betti(small_worlds1[k], [0.1], 1) for k in range(len(small_worlds1))]
betti_ns1=np.array([[betti_ns1[i][j][0] for j in range(len(betti_ns1[0]))]for i in range(len(betti_ns1))])
betti_mean=np.mean(betti_ns1,axis=0)
probs=np.arange(100)*0.01
betti_sorted=np.sort(betti_ns1,axis=0)
#take quantile
lower=betti_sorted[int(np.floor(5*np.shape(betti_sorted)[0]/100))]
upper=betti_sorted[int(np.floor(95*np.shape(betti_sorted)[0]/100))]
fig,ax=plt.subplots()
ax.plot(np.arange(len(betti_mean))*0.01,betti_mean)
ax.fill_between(np.arange(len(betti_mean))*0.01,lower,upper,alpha=0.2)
ax.set_xlabel('p')
ax.set_ylabel('1-betti number')


small_worlds2=[[generate_small_world_correlation_matrix(90,0.001*i, k) for i in range(200)] for k in range(100)]

betti_ns=[betti(small_worlds2[k], [0.1], 1) for k in range(len(small_worlds2))]
betti_ns=np.array([[betti_ns[i][j][0] for j in range(len(betti_ns[0]))]for i in range(len(betti_ns))])
betti_mean=np.mean(betti_ns,axis=0)
probs=np.arange(200)*0.001
betti_sorted=np.sort(betti_ns,axis=0)
#take quantile
lower=betti_sorted[int(np.floor(5*np.shape(betti_sorted)[0]/100))]
upper=betti_sorted[int(np.floor(95*np.shape(betti_sorted)[0]/100))]
fig,ax=plt.subplots()
ax.plot(probs,betti_mean)
ax.fill_between(probs,lower,upper,alpha=0.2)
ax.set_xlabel('p')
ax.set_ylabel('1-Betti number')


#ora ho stampato i grafici dei betti
# quello che ha senso vedere ora sono le proprietà non topologiche di questi network
from bct import efficiency_bin

ring=generate_ring_lattice_correlation_matrix(90)
random=generate_Erdos_Renyi_correlation_matrix(90,0)
WS=generate_small_world_correlation_matrix(90, 0.015)
cm_brain = np.loadtxt('regional-differentiation-based-on-graph-nodal-statistics-for-functional-brain-connectivity-networks-characterization/DATA/cor_mat_coma/cor_mat_90_wave_scale_3_ts_Control_grey_matter_ROI_90_ts_Control_0.txt.txt',usecols=range(1,91), skiprows=1, dtype='float')

adj_matrices_ring=[adj_matrix_notconncted(ring, eps) for eps in np.linspace(0, 1,100)]
glob_eff_ring=[efficiency_bin(matrix) for matrix in adj_matrices_ring]
adj_matrices_random=[adj_matrix_notconncted(random, eps) for eps in np.linspace(0, 1,100)]
glob_eff_random=[efficiency_bin(matrix) for matrix in adj_matrices_random]
adj_matrices_WS=[adj_matrix_notconncted(WS, eps) for eps in np.linspace(0, 1,100)]
glob_eff_WS=[efficiency_bin(matrix) for matrix in adj_matrices_WS]
adj_matrices_brain=[adj_matrix_notconncted(cm_brain, eps) for eps in np.linspace(0, 1,100)]
glob_eff_brain=[efficiency_bin(matrix) for matrix in adj_matrices_brain]

fig,ax=plt.subplots()
ax.plot(np.linspace(0, 1,100),glob_eff_ring)
ax.plot(np.linspace(0, 1,100),glob_eff_random)
ax.plot(np.linspace(0, 1,100),glob_eff_WS)
ax.plot(np.linspace(0, 1,100),glob_eff_brain)
loc_eff_ring=[np.mean(efficiency_bin(matrix,1)) for matrix in adj_matrices_ring]
loc_eff_random=[np.mean(efficiency_bin(matrix,1)) for matrix in adj_matrices_random]
loc_eff_WS=[np.mean(efficiency_bin(matrix,1)) for matrix in adj_matrices_WS]
loc_eff_brain=[np.mean(efficiency_bin(matrix,1)) for matrix in adj_matrices_brain]
fig,ax=plt.subplots()
ax.plot(np.linspace(0, 1,100),loc_eff_ring)
ax.plot(np.linspace(0, 1,100),loc_eff_random)
ax.plot(np.linspace(0, 1,100),loc_eff_WS)
ax.plot(np.linspace(0, 1,100),loc_eff_brain)


np.mean(efficiency_bin(network,1))
nx.draw(nx.from_numpy_matrix(network))


###Per ora abbiamo disegnato solo lale efficenze locali e globali dei cervelli.
# Vediamo se è possibile fare inferenza su un paramentro p sui cervelli
# dobbiamo costrire uno stimatore. La cosa più semplice è prendere punti nello
# small worls regime (0.05-0.34) e vedere i parametri di un polinomio calcolato 
# come variano in funzione di p. In pratica è abbiamo due attributi su cui fare 
# regressione polinomiale: costo e probabilità
# sta cosa è una stronzata: la probabilità è lattributo su cui fare regressione
# i valori dei vari costi sono gli attributi indipendenti

from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.metrics import mean_squared_error

small_worlds1=[[generate_small_world_correlation_matrix(90,0.05*i, k) for i in range(20)] for k in range(10)]
small_worlds1=np.array(small_worlds1)
ws_corr=np.reshape(small_worlds1,[200,90,90])
probabilites=np.array([[0.05*i for i in range(20)] for k in range(10)])
probs=probabilites.flatten()
adj_matrices=np.array([[adj_matrix_notconncted(matrix, eps) for eps in np.linspace(0.05, 0.34,10)] for matrix in ws_corr])
glob_eff=np.array([[efficiency_bin(matrix) for matrix in single_ws] for single_ws in adj_matrices])
X=np.concatenate((glob_eff,glob_eff**2,glob_eff**3),axis=1)

np.corrcoef(X.T)


X_train, X_test, y_train, y_test = train_test_split(X, probs, test_size=0.2, random_state=0)
reg=linear_model.LinearRegression()

reg.fit(X_train,y_train)
reg.coef_
y_pred=reg.predict(X_test)
mean_squared_error(y_test,y_pred)
ind=np.argsort(y_test)
plt.plot(y_test[ind],y_pred[ind])


adj_matrices_brain=[adj_matrix_notconncted(cm_brain, eps) for eps in np.linspace(0.05, 0.34,10)]
glob_eff_brain=np.array([efficiency_bin(matrix) for matrix in adj_matrices_brain])
reg.predict([np.concatenate((glob_eff_brain,glob_eff_brain**2,glob_eff_brain**3))])

# non mi piace molto come viene qua, quello che si potrebbe fare è: approssimare ogni curva a un polinomio e mettere come regressori quei coefficienti
adj_matrices_brain=[adj_matrix_notconncted(cm_brain, eps) for eps in np.linspace(0.05, 0.34,10)]
glob_eff_brain=np.array([efficiency_bin(matrix) for matrix in adj_matrices_brain])
costs=np.reshape(np.linspace(0.05, 0.34,10),[10,1])
X=np.concatenate((costs,costs**2,costs**3),axis=1)
reg.fit(X,glob_eff_brain)
a,b,c=reg.coef_
d=reg.intercept_
x=np.linspace(0.05, 0.34,100)
y=a*x+b*x**2+c*x**3+d
fig,ax=plt.subplots()
ax.plot(x,y)
ax.plot(np.linspace(0.05, 0.34,10),glob_eff_brain)

#adesso che abbiamo capito come si fa, proviamo a farlo su ogni cervello e a fare regressione sui parametri a,b,c,d
reg_interna=linear_model.LinearRegression(fit_intercept=False)
costs=np.reshape(np.linspace(0.05, 0.34,10),[10,1])
X_interna=np.concatenate((costs,costs**2,costs**3),axis=1)
X_esterna=[]
for corr_matrix in ws_corr:
    #faccio la regressione per i valori a,b,c,d
    adj_matrices=[adj_matrix_notconncted(corr_matrix, eps[0]) for eps in costs]
    glob_eff=np.array([efficiency_bin(matrix) for matrix in adj_matrices])
    reg_interna.fit(X_interna,glob_eff)
    #X_esterna.append([reg_interna.intercept_,reg_interna.coef_[0],reg_interna.coef_[1],reg_interna.coef_[2]])
    X_esterna.append([reg_interna.coef_[0],reg_interna.coef_[1],reg_interna.coef_[2]])

X_esterna=np.array(X_esterna)    
reg_esterna=linear_model.LinearRegression(fit_intercept=False)   
X_train, X_test, y_train, y_test = train_test_split(X_esterna, probs, test_size=0.2, random_state=0)
reg_esterna.fit(X_train,y_train)
reg_esterna.coef_
y_pred=reg_esterna.predict(X_test)
mean_squared_error(y_test,y_pred)
ind=np.argsort(y_test)
plt.plot(y_test[ind],y_pred[ind])

adj_matrices_brain=[adj_matrix_notconncted(cm_brain, eps) for eps in np.linspace(0.05, 0.34,10)]
glob_eff_brain=np.array([efficiency_bin(matrix) for matrix in adj_matrices_brain])
costs=np.reshape(np.linspace(0.05, 0.34,10),[10,1])
X=np.concatenate((costs,costs**2,costs**3),axis=1)
reg_interna.fit(X,glob_eff_brain)
a,b,c=reg_interna.coef_
d=reg_interna.intercept_

#adj_matrices_brain=[adj_matrix_notconncted(cm_brain, eps) for eps in np.linspace(0.05, 0.34,10)]
#glob_eff_brain=np.array([efficiency_bin(matrix) for matrix in adj_matrices_brain])
reg_esterna.predict([[a,b,c]])

#sta roba è venuta un po' una schifezza non mi dà risultati buoni. Bisogna ragionare un bel po'






#proviamo a calcolare la vietoris risps homology dato distanze sul network. Prendo un esempio di WS
network=small_worlds2[0][15]
network=adj_matrix_notconncted(network, 0.1)
G=nx.from_numpy_matrix(network)
distances=dict(nx.shortest_path_length(G))
dist=np.array([[distances[j][i] for i in range(90)]for j in range(90)])


from gtda.homology import VietorisRipsPersistence

VR = VietorisRipsPersistence(metric="precomputed",homology_dimensions=[1,2])
diagrams = VR.fit_transform([dist])

#laconfronto con quello fatto da noi
homo=homology(small_worlds2[0][15],[0.1],1)

from gtda.plotting import plot_diagram
plot_diagram(diagrams[0])
#non mi plotta le cose non ne capisco il motivo



#ora mapper
import numpy as np
import pandas as pd  # Not a requirement of giotto-tda, but is compatible with the gtda.mapper module

# Data viz
from gtda.plotting import plot_point_cloud

# TDA magic
from gtda.mapper import (
    CubicalCover,
    make_mapper_pipeline,
    Projection,
    plot_static_mapper_graph,
    plot_interactive_mapper_graph,
    MapperInteractivePlotter
)

# ML tools
from sklearn import datasets
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA

data, _ = datasets.make_circles(n_samples=5000, noise=0.05, factor=0.3, random_state=42)

plot_point_cloud(data)
filter_func = Projection(columns=[0, 1])
cover = CubicalCover(n_intervals=10, overlap_frac=0.3)
clusterer = DBSCAN()
n_jobs = 1
pipe = make_mapper_pipeline(
    filter_func=filter_func,
    cover=cover,
    clusterer=clusterer,
    verbose=False,
    n_jobs=n_jobs,
)
fig = plot_static_mapper_graph(pipe, data)
fig.show(config={'scrollZoom': True})



'''

