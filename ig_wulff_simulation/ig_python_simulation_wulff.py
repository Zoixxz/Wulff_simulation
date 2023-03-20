import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import ig_percoSimAux as aux
import ig_percoSimIteration as upd
import math
import time as t
print("Python version")
print(sys.version)
print("Version info.")
print(sys.version_info)
print("igraph version: " + ig.__version__)


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


N = 40  # we build a box from -N to N
# m is size of inside box, must be strictly smaller than N
m = 10
p = 0.7


iterations = 300001
step_size = 10  # step size for boundary size plotting
draw_step_size = 50000  # step size of drawing
time_calc = 1000

folderName = ('G_figures/' + 'non_rev_G_simulation_N' + str(N) + 'm' + str(m)
              + 'iterations' + str(iterations) + 'p' + str(p) + 'EXP')
createFolder(folderName)
filetype = '.png'  # use the filetype used to store the plots of the graph
print('iterations: ', iterations)


"""
G = Graph as built by gt
"""

G = ig.Graph((2 * N + 1) ** 2)
G.simplify()  # initialise the ig graph with vertices

# V contains the layout of the graph for plotting
G, V, interior_edges, boundary_edges = aux.createConfigNorm(
    G, N, m, p)
total_edges = len(G.es)
print('Total number of edges in graph = ', total_edges)
# ig.plot(G, layout=V, target='myfile.pdf')

boundary_size = []
start = t.time()
for i in range(iterations):

    G, boundary_edges, interior_edges = upd.wulffIteration_edges(
        G, boundary_edges, interior_edges, p, N)
    if(i % 10 == 0):
        print('iteration number ', i)
        if(i == time_calc):
            end = t.time()
            print('total time needed [s], ', (iterations / time_calc)
                  * (end - start))

    if(i % step_size == 0):
        boundary_size.append(len(boundary_edges))

    if(i % draw_step_size == 0):
        fig, ax = plt.subplots(figsize=(10, 10))
        ig.plot(G, layout=V, target=ax, vertex_size=0)
        plt.title('iteration number ' + str(i)
                  + ', |C| = ' + str(total_edges)
                  + ', p = ' + str(p), fontsize=20)
        plt.savefig(folderName + '/' + str(i) + filetype)
        plt.clf()
        plt.close()

plt.plot(np.arange(1, iterations + 1, step_size), boundary_size[:])
plt.xlabel('steps')
plt.ylabel('boundary size')
plt.savefig(folderName + '/' + 'boundary_size' + '_iterations_'
            + str(iterations) + 'EXP' + '.pdf')
plt.clf()
plt.close()
