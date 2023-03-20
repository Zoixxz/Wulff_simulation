from scipy.stats import bernoulli
import percoSimIterationAux as aux


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
    empty array of Vertices
    build the vertex set [a list of vertices]
    V is a 2N+1 time 2N+1 array
    V[0] = [-N, -N] and V[1] = [-N+1, -N]
    """
    V = []

    for j in range(- N, N + 1):
        for i in range(- N, N + 1):
            V.append([i, j])  # add the vertices of the graph
    V = [tuple(v) for v in V]
    # add all the vertices/nodes to G
    for v in V:
        G.add_node(tuple(v), pos=tuple(v))  # convert to tuple to be iterable

    # add all the eges (random) to G
    # iterate over horizontal (x-dir) edges
    # add random background edges to the graph
    for i in range((2 * N + 1) * (2 * N + 1) - 1):
        if(isNeighbour(V[i], V[i + 1]) and bernoulli.rvs(p) == 1):
            G.add_edge(V[i], V[i + 1])

    # iterate over vertical (y-dir) edges
    # add random background edges to the graph
    for i in range((2 * N + 1) * 2 * N):
        if(isNeighbour(V[i], V[i + (2 * N + 1)]) and bernoulli.rvs(p) == 1):
            G.add_edge(V[i], V[i + (2 * N + 1)])

    # add the edges for initial configuration
    # initialise the two edge sets
    interior_edges = []  # set of sets of edges
    boundary_edges = []  # set of sets of edges
    for v in V:
        neighbours_v = aux.neighbours(v)
        for w in neighbours_v:
            if(graphDistanceOne(v) <= m):  # takes care of bdry automatically
                if(graphDistanceOne(w) <= m):  # interior edge
                    G.add_edge(v, w)
                    interior_edges.append(tuple([v, w]))
                if(graphDistanceOne(w) == m + 1):  # boundary edge
                    if(G.has_edge(v, w)):
                        G.remove_edge(v, w)
                    boundary_edges.append(tuple([v, w]))
    return G, V, interior_edges, boundary_edges
