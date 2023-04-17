import os
import numpy as np
import matplotlib.pyplot as plt
import igraph as ig
import ig_percoSimAux as aux
import ig_percoSimIterationAux as itaux
import ig_percoSimIteration as upd

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


def performSimulation(v):
    """
    INPUT
    N               int         parameter of box size (half side length +1)
    m               int         less than N, distance for norm for initial
                                configuration
    p               float       number in (0, 1), perco density
    iterations      int         number of iterations which are to be
                                performed
    draw_step_size  int         number of steps after which the graph is
                                drawn (more pictures cost more time)
    do_plots        bool        bolean if True pictures will be plotted

    PROCEDURE
    This function performs a simulation with the above specified paremeters
    and creates a folder in the cwd where it stores the plotted figures
    graphs are stored as a png file and plots as a pdf.
    The filetype for the plots of the clusters can be changes via the filetype
    string.

    RETURNS
    time_series_boundary        return the time series of the boundary
                                size for each simulation step
    """
    N = v[0]  # we build a box from -N to N
    m = v[1]  # m is size of inside box, must be strictly smaller than N
    p = v[2]
    iterations = v[3]

    draw_step_size = v[4]  # step size of drawing
    do_plots = v[5]
    j = v[6]
    if(len(v) != 7):
        return print('ERROR: wrong size of input value for performSimulation')

    folderName = ('G_figures/' + 'rev_G_simulation_N' + str(N) + 'm' + str(m)
                  + 'iterations' + str(iterations) + 'p' + str(p) + '/'
                  + 'number' + str(j))
    createFolder(folderName)
    filetype = '.png'  # use the filetype used to store the plots of the graph
    # print('iterations: ', iterations)

    """
    initialise the lists into which we store the values of the iteration
    """
    N_C_list = np.zeros(iterations)
    N_D_list = np.zeros(iterations)
    TP_list = np.zeros(iterations)
    TP_non_list = np.zeros(iterations)
    time_series_boundary = np.zeros(iterations)

    """
    G = Graph as built by ig
    """
    G = ig.Graph((2 * N + 1) ** 2)  # initialise the ig graph with vertices
    G.simplify()  # we want a simple graph (not directed)

    # V contains the layout of the graph for plotting (vertex coordinates)
    G, V, interior_edges, boundary_edges = aux.createConfigNorm(
        G, N, m)
    total_edges = len(G.es)
    # print('Total number of edges in graph = ', total_edges)
    # print('boundary size = ', len(boundary_edges))

    for t in range(iterations):

        G, boundary_edges, interior_edges, N_C, N_D, TP, TP_non = upd.wulffIteration_edges(
            G, boundary_edges, interior_edges, p, N)
        N_C_list[t] = N_C
        N_D_list[t] = N_D
        TP_list[t] = TP
        TP_non_list[t] = TP_non

        time_series_boundary[t] = len(boundary_edges)

        if(t % draw_step_size == 0 and do_plots):
            fig, ax = plt.subplots(figsize=(10, 10))
            F = G.copy()
            for e_frozenset in boundary_edges:
                e = [c for c in e_frozenset]
                F.add_edge(itaux.coord_to_array(
                    N, e[0]), itaux.coord_to_array(N, e[1]))

            ig.plot(F, layout=V, target=ax, vertex_size=0,
                    edge_color='red', width=0.02)
            ig.plot(G, layout=V, target=ax, vertex_size=0)

            plt.title('iteration number ' + str(t)
                      + ', |C| = ' + str(total_edges)
                      + ', p = ' + str(p), fontsize=20)
            plt.savefig(folderName + '/' + str(t) + filetype)
            plt.clf()
            plt.close()

    """
    Create the plots comparing N(C) and N(D)
    first plot: compares the N(C) and N(D)
    second plot: compares the transition probabilities to the non reversible
    proposed chain
    """

    s = 1  # s specifies the size of the dots in the plots (s = 1) works well
    fig, ax = plt.subplots(2)
    x_axis = np.arange(len(N_C_list))
    ax[0].scatter(x_axis, N_C_list, s=s, color='blue')
    ax[0].scatter(x_axis, N_D_list, s=s, color='red')
    ax[0].set_title(r'plot of $N(C)$ and $N(D)$')
    ax[0].set_ylabel(r'$\frac{N(C)}{N(D)}$')
    ax[0].legend((r'$N(C)$', r'$N(D)$'))
    quotient = [x / y for (x, y) in zip(N_C_list, N_D_list)]
    x_axis = np.arange(len(quotient))
    ax[1].scatter(x_axis, quotient, s=s, color='blue')
    ax[1].set_title(r'plot of $\frac{N(C)}{N(D)}$')
    ax[1].set_ylabel(r'$\frac{N(C)}{N(D)}$')
    fig.tight_layout()
    plt.savefig(folderName + '/' + 'NCNDplot' + '.pdf')
    plt.clf()
    plt.close()

    fig, ax = plt.subplots(2)
    x_axis = np.arange(len(TP_list))
    ax[0].scatter(x_axis, TP_list, s=s, color='blue')
    ax[0].scatter(x_axis, TP_non_list, s=s, color='red')
    ax[0].set_title('plot of transitioning probabilities: values')
    ax[0].set_ylabel('TPs')
    ax[0].legend(('rev', 'non-rev'))
    difference = [x - y for (x, y) in zip(TP_list, TP_non_list)]
    x_axis = np.arange(len(difference))
    ax[1].scatter(x_axis, difference, s=s, color='blue')
    ax[1].set_title('TPs: difference (rev - non-rev)')
    ax[1].set_ylabel('difference TPs')
    fig.tight_layout()
    plt.savefig(folderName + '/' + 'TP_plots' + '.pdf')
    plt.clf()
    plt.close()

    """
    Here we create the boundary size plots over the entire duration of the
    simulation with time steps of length draw_step_size
    """
    x_axis = np.arange(len(time_series_boundary))
    plt.scatter(x_axis, time_series_boundary, s=s)
    plt.xlabel('steps')
    plt.ylabel('boundary size')
    plt.savefig(folderName + '/' + 'boundary_size' + '_iterations_'
                + str(iterations) + '.pdf')
    plt.clf()
    plt.close()

    return time_series_boundary
