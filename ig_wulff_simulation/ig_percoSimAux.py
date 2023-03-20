from scipy.stats import bernoulli
import igraph as ig
import ig_percoSimIterationAux as itaux



"""
This file contains auxiliary functions which are needed to create
the initial configuration

Functions:
     isNeighbour
        returns whether two vectors/coordinates are neighbors in
        the lattice Z^2

     graphDistanceMax(v)
        returns the maximum-norm of a vector/coordinate

     relativePosition Boundary
        returns a string containing the info which part of the boundary
        a vertex/coordinate is in

    createConfigSquare
"""


def isNeighbour(A, B):
    # return true/falso is A is neighbour of B
    # for A and B arrays with length 2
    if(len(A) == len(B)):
        diff = [(A[i] - B[i]) for i in range(2)]
        distance = abs(sum(diff))
        if(A != B and distance == 1):
            return 1
        else:
            return 0
    else:
        print('lengths do not match')


def graphDistanceMax(v):
    # v: vertex with coordinates
    # returns max norm of the vertes/coordinate
    # works is any dimension
    dist = max([abs(x) for x in v])
    return dist


def graphDistanceOne(v):
    dist = sum([abs(x) for x in v])
    return dist


def createConfigNorm(G, N, m, p):
    """
    G an ig graph with (2 * N + 1) ** 2 vertices
    empty array of Vertices
    build the vertex set [a list of vertices]
    V is a 2N+1 time 2N+1 array
    V[0] = [-N, -N] and V[1] = [-N, -N]
    """

    V = []
    for j in range(- N, N + 1):  # y-coordinate index of vertex
        for i in range(- N, N + 1):  # x-coordinated index of vertex
            V.append((i, j))  # add the vertices of the graph

    # add the edges for initial configuration
    # initialise the interior and boundary edge sets
    interior_edges = []  # set of sets of edges (tuples with coordinates)
    boundary_edges = []  # set of sets of edges (tuples with coordinates)
    for v in V:
        neighbours_v = itaux.neighbours(v)
        for w in neighbours_v:
            if(graphDistanceOne(v) <= m):  # takes care of bdry autom.
                if(graphDistanceOne(w) <= m):  # interior edge
                    if(not G.are_connected(itaux.coord_to_array(N, v),
                                           itaux.coord_to_array(N, w))):
                        e = G.add_edge(itaux.coord_to_array(N, v),
                                       itaux.coord_to_array(N, w))
                    interior_edges.append((v, w))
                if(graphDistanceOne(w) == m + 1):  # boundary edge
                    if(G.are_connected(itaux.coord_to_array(N, v),
                                       itaux.coord_to_array(N, w))):
                        e = G.get_eid(itaux.coord_to_array(N, v),
                                      itaux.coord_to_array(N, w))
                        G.delete_edges([e])
                    boundary_edges.append((v, w))
    return G, V, interior_edges, boundary_edges
