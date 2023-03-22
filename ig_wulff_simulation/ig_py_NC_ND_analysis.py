import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import ig_percoSimAux as aux
import ig_percoSimIteration as upd
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


N = 20  # we build a box from -N to N
# m is size of inside box, must be strictly smaller than N
m = 10
p = 0.55


iterations = 1000
step_size = 100  # step size for boundary size plotting
draw_step_size = 400000  # step size of drawing
time_calc = 1000

folderName = ('G_figures/' + 'N(C)_N(D)_comparison')
createFolder(folderName)
# use the filetype used to store the plots of the graph
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

N_C_list = []
N_D_list = []


for i in range(iterations):
    if(i % 10 == 0):
        print('iteration ', i)

    G, boundary_edges, interior_edges, N_C, N_D = upd.wulffIteration_edges(
        G, boundary_edges, interior_edges, p, N)
    N_C_list.append(N_C)
    N_D_list.append(N_D)
x_axis = np.arange(len(N_C_list))
plt.plot(x_axis, N_C_list, '.', label='N(C)')
plt.plot(x_axis, N_D_list, '.', label='N(D)')
plt.xlabel('steps')
plt.ylabel('N_size')
plt.title('iter=' + str(iterations) + ' m = ' + str(m) + ' N = ' + str(N))
plt.legend()
plt.savefig(folderName + '/' + 'plot_neighbours.pdf')
plt.clf()
plt.close()

quotient = [x / y for (x, y) in zip(N_C_list, N_D_list) if x != 0 and y != 0]
x_axis = np.arange(len(quotient))
plt.plot(x_axis, quotient, '.', label='quotient')
plt.xlabel('steps')
plt.ylabel('quotient')
plt.title('plot of quotient')
plt.legend()
plt.savefig(folderName + '/' + 'plot_neighbours_quotient.pdf')
plt.clf()
plt.close()

