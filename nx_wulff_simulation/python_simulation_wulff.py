import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import percoSimAux as aux
import percoSimIteration as upd
import math
import sys
# import time as t
print("Python version")
print(sys.version)
print("Version info.")
print(sys.version_info)
print("networkx version: " + nx.__version__)


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


"""
set simulation parameters
N = side length of box which will be simulated
m = cluster size (must be m <= N)
p = density parameter
"""


def orderOfMagnitude(number):
    return math.floor(math.log(number, 10))


N = 30  # we build a box from -N to N
# m is size of inside box, must be strictly smaller than N
m = 15
p = 0.6


iterations = 100001
step_size = 10
draw_step_size = 1000

folderName = ('G_figures/' + '_G_simulation_N' + str(N) + 'm' + str(m)
              + 'iterations' + str(iterations))
createFolder(folderName)
filetype = '.png'  # use the filetype used to store the plots of the graph
print('iterations: ', iterations)


"""
G = Graph as built by nx
"""

# initialize G as the graph
G = nx.Graph()
G, V, interior_edges, boundary_edges = aux.createConfigNorm(G, N, m, p)
boundary_size = []
for i in range(iterations):
    G, boundary_edges, interior_edges, bridge_edges = upd.wulffIteration_edges(
        G, boundary_edges, interior_edges, p)

    if(i % 10 == 0):
        print('iteration number ', i)

    if(i % step_size == 0):
        boundary_size.append(len(boundary_edges))

    if(i % draw_step_size == 0):

        Cluster = G.edge_subgraph(tuple(interior_edges))
        boundary_subgraph = G.edge_subgraph(tuple(boundary_edges))
        # bridges_graph = G.edge_subgraph(tuple(bridge_edges))
        bdry_subgr = boundary_subgraph.copy()

        for v in boundary_edges:
            bdry_subgr.add_edge(v[0], v[1])

        fig = plt.figure(figsize=(10, 10))

        plt.title(str(i), fontsize=20)

        posG = nx.get_node_attributes(G, 'pos')
        posC = nx.get_node_attributes(Cluster, 'pos')
        posBdry = nx.get_node_attributes(bdry_subgr, 'pos')
        # pos_bridges = nx.get_node_attributes(bridges_graph, 'pos')

        nx.draw(G, pos=posG, node_size=(1 / N), width=2)
        nx.draw(Cluster, pos=posC, node_size=(1 / N),
                width=3, edge_color='blue')
        nx.draw(bdry_subgr, pos=posBdry, node_size=(1 / N),
                edge_color='red', width=1)
        # nx.draw(bridges_graph, pos=pos_bridges, node_size=(1/N), edge_color=
        # 'green', width=7)
        
        plt.savefig(folderName + '/' + str(i) + filetype)
        plt.clf()
        plt.close()

plt.plot(np.arange(1, iterations + 1, step_size), boundary_size[:])
plt.xlabel('steps')
plt.ylabel('boundary size')
plt.savefig(folderName + '/' + 'boundary_size' + '_iterations_'
            + str(iterations) + '.pdf')
plt.clf()
plt.close()
