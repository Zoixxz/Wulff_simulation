import os
import numpy as np
import igraph as ig
import ig_percoSimAux_runtime as aux
import ig_percoSimIteration_runtime as upd
import time as ti

"""
COMMENTED OUT>>> left in for potential debugging purposes
the following print statements give info on the versions of
python used and the igraph library
print("Python version")
print(sys.version)
print("Version info.")
print(sys.version_info)
print("igraph version: " + ig.__version__)

"""


def createFolder(directory):
    """
    creates a folder in the directory path
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def performSimulation(input_values):
    """
    INPUT
    N               int         parameter of box size (half side length +1)
    m               int         less than N, distance for norm for initial
                                configuration
    p               float       number in (0, 1), perco density
    iterations      int         number of iterations which are to be
                                performed

    PROCEDURE
    This function performs a simulation with the above specified paremeters
    and creates a folder in the cwd where it stores the plotted figures
    graphs are stored as a png file and plots as a pdf.
    The filetype for the plots of the clusters can be changes via the filetype
    string.

    RETURNS
    time_measurements           time series of the measurements of each
                                iteration step
    """
    if(len(input_values) != 4):
        return print('ERROR: wrong size of input value for performSimulation')
    N = input_values[0]
    m = input_values[1]
    p = input_values[2]
    iterations = input_values[3]
    # initialise the lists into which we store the values of the iteration

    time_series_boundary = np.zeros(iterations)
    time_measurements = np.zeros(iterations)

    """
    G = Graph as built by ig
    """
    G = ig.Graph((2 * N + 1) ** 2)  # initialise the ig graph with vertices
    G.simplify()  # we want a simple graph (not directed)

    # V contains the layout of the graph for plotting (vertex coordinates)
    G, V, interior_edges, boundary_edges = aux.createConfigNorm(
        G, N, m)
    for t in range(iterations):
        start = ti.time()
        G, boundary_edges, interior_edges, N_C, N_D, TP, TP_non = upd.wulffIteration_edges(
            G, boundary_edges, interior_edges, p, N)
        end = ti.time()

        time_step = end - start
        time_measurements[t] = time_step

        time_series_boundary[t] = len(boundary_edges)

    return time_measurements
