import ig_percoSimIterationAux as itaux


"""
This file contains auxiliary functions which are needed to create
the initial configuration

Functions:


    graphDistanceMax(v)
        returns the maximum-norm of a vector/coordinate

    graphDistanceOne(v)
        returns the one-norm of a vector/coordinate

    createConfigNorm(G, N, m, p)
        creates the initial configuration using one of the two norms above
        returns Graph, Vertex coordinate list, list of interior edges, list of
                boundary edges
"""


def graphDistanceMax(v):
    """
    INPUT:
        v      tuple of ints       vertex with coordinates
    OUTPUT:
        dist   int                 max-norm of the vertes/coordinate
    """
    dist = max([abs(x) for x in v])
    return dist


def graphDistanceOne(v):
    """
    INPUT:
        v      tuple of ints       vertex with coordinates
    OUTPUT:
        dist   int                 one-norm of the vertes/coordinate
    """
    dist = sum([abs(x) for x in v])
    return dist


def createConfigNorm(G, N, m):
    """
    INPUT:
        G               ig graph            with (2 * N + 1) ** 2 vertices
        N               int                 half side length of box
        m               int                 distance for norm for config

    build the vertex set [a list of vertices]
    V is a 2N+1 time 2N+1 array
    V[0] = [-N, N] and V[1] = [-N, N]

    OUTPUT:
        G               ig graph            Graph initialised with edges
        V               list of tuples      list with coordinates for plotting
        interior_edges  list of tuples      list with edges of interior
        boundary_edges  list of tuples      list with edges of boundary
    """

    V = []
    for j in range(- N, N + 1):  # y-coordinate index of vertex
        for i in range(- N, N + 1):  # x-coordinated index of vertex
            V.append((i, j))  # add the vertices of the graph

    # add the edges for initial configuration
    # initialise the interior and boundary edge sets
    interior_edges = []  # set of sets of edges (frozensets with coordinates)
    boundary_edges = []  # set of sets of edges (frozensets with coordinates)
    for v in V:
        neighbours_v = itaux.neighbours(v)
        for w in neighbours_v:
            if(graphDistanceOne(v) <= m):  # takes care of bdry autom.
                if(graphDistanceOne(w) <= m):  # interior edge
                    if(not G.are_connected(itaux.coord_to_array(N, v),
                                           itaux.coord_to_array(N, w))):
                        e = G.add_edge(itaux.coord_to_array(N, v),
                                       itaux.coord_to_array(N, w))
                    interior_edges.append(frozenset((v, w)))
                if(graphDistanceOne(w) == m + 1):  # boundary edge
                    if(G.are_connected(itaux.coord_to_array(N, v),
                                       itaux.coord_to_array(N, w))):
                        e = G.get_eid(itaux.coord_to_array(N, v),
                                      itaux.coord_to_array(N, w))
                        G.delete_edges([e])
                    boundary_edges.append(frozenset((v, w)))
    return G, V, interior_edges, boundary_edges
