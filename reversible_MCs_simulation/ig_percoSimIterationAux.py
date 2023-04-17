def neighbours(x):
    """
    takes x [tuple] a vertex as input and return
    a list of all neighbours in Z^2
    """
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
            for i in range(len(iterate_neighbours))]


def neighbours_vid(N, y):
    """
    takes y [tuple] a vertex as input and return
    a list of all neighbours (in Z^2 as a grpah) with the vertex coordinates
    """
    x = array_to_coord(N, y)
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    V = [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
         for i in range(len(iterate_neighbours))]
    return [coord_to_array(N, v) for v in V]


def transition_rate_edges(p, Delta_boundary):
    """
    computes a transition rate
    """
    return (1 - p) ** Delta_boundary


def transition_proba_rev(G, p, N, boundary_edges,
                         boundary_edges_updated,
                         e_bdry, e_interior):
    """
    INPUT:
        G                       ig Graph
        p                       density parameter bond percolation
        interior_edges          old interior edges of the cluster
        boundary_edges          old boundary edges of the cluster
        interior_edges_updated  new interior edges of the cluster
        boundary_edges_updated  new boundary edges of the cluster
        e_bdry                  boundary edge
        e_interior              interior edge

    function computes the reversible transition probability to go
    from G to to the updated G by exchanging the boundary edges
    and the interior edge

    returns:
        TP                      probability to transition to the other state
    """
    edge_count = G.ecount()
    F = G.copy()
    N_C = 0
    N_D = 0

    bridges_no_change = F.es[F.bridges()]
    bridges_number_no_change = len([e for e in bridges_no_change
                                    if not isLeaf_edge_id(F, e)])
    # compute N(C), the neighbouring states of C for the relation R
    for w in boundary_edges:  # w is a frozenset
        v = [c for c in w]  # unpack frozenset to make iterable
        F.add_edge(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        e_id = F.get_eid(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        e_es = F.es[e_id]
        # check if leaf
        if(isLeaf_edge(N, F, v)):
            v1, v2 = e_es.source, e_es.target
            if(F.degree(v1) >= 3 or F.degree(v2) >= 3):
                # edge which was added is not adjacent to leaf in G
                N_C += edge_count - bridges_number_no_change
                e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                      coord_to_array(N, v[1]))
                F.delete_edges([e_bdry_id])
                continue
            if(F.degree(v1) <= 2 and F.degree(v2) <= 2):
                # edge which was added is adjacent to leaf in G
                N_C += edge_count - (bridges_number_no_change + 1)
                e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                      coord_to_array(N, v[1]))
                F.delete_edges([e_bdry_id])
                continue
        else:
            bridges = F.es[F.bridges()]
            bridges_number = len([e for e in bridges
                                  if not isLeaf_edge_id(F, e)])
            N_C += edge_count - bridges_number
            e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                  coord_to_array(N, v[1]))
            F.delete_edges([e_bdry_id])
            continue

    # compute N(D), the neighbouring states of D for the relation R
    # update the graph F so that we have configuaration D
    F.add_edge(coord_to_array(N, e_bdry[0]), coord_to_array(N, e_bdry[1]))
    e_id = F.get_eid(coord_to_array(N, e_interior[0]),
                     coord_to_array(N, e_interior[1]))
    F.delete_edges([e_id])

    bridges_no_change_D = F.es[F.bridges()]
    bridges_number_no_change_D = len([e for e in bridges_no_change_D
                                      if not isLeaf_edge_id(F, e)])
    # compute N(C), the neighbouring states of C for the relation R
    for w in boundary_edges:  # w is a frozenset
        v = [c for c in w]
        F.add_edge(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        e_id = F.get_eid(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        e_es = F.es[e_id]
        # check if leaf
        if(isLeaf_edge(N, F, v)):
            v1, v2 = e_es.source, e_es.target
            if(F.degree(v1) >= 3 or F.degree(v2) >= 3):
                # edge which was added is not adjacent to leaf in G
                N_D += edge_count - bridges_number_no_change_D
                e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                      coord_to_array(N, v[1]))
                F.delete_edges([e_bdry_id])
                continue
            if(F.degree(v1) <= 2 and F.degree(v2) <= 2):
                # edge which was added is adjacent to leaf in G
                N_D += edge_count - (bridges_number_no_change_D + 1)
                e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                      coord_to_array(N, v[1]))
                F.delete_edges([e_bdry_id])
                continue
        else:
            bridges = F.es[F.bridges()]
            bridges_number = len([e for e in bridges
                                  if not isLeaf_edge_id(F, e)])
            N_D += edge_count - bridges_number
            e_bdry_id = F.get_eid(coord_to_array(N, v[0]),
                                  coord_to_array(N, v[1]))
            F.delete_edges([e_bdry_id])
            continue

    Delta_boundary = len(boundary_edges_updated) - len(boundary_edges)
    """
    return ((((1 - p) ** Delta_boundary) * (N_C / N_D)) / (1 + ((1 - p) ** Delta_boundary) * (N_C / N_D))), N_C, N_D
    """
    return min((1 - p) ** Delta_boundary * (N_C / N_D), 1), N_C, N_D


def array_to_coord(N, location):
    """
    function computes the coordinate of vertex in G given location in array

    INPUT
        N        int     size parameter of box
        location int     position in array

    OUTPUT
        coord    tuple   coordinate of vertex in graph
    """
    m_y, m_x = divmod(location, 2 * N + 1)
    return (m_x - N, m_y - N)


def coord_to_array(N, coordinate):
    """
    function computes the coordinate of vertex in G given location in array

    INPUT
        N        int     size parameter of box
        location tuple   coordinate of vertex

    OUTPUT
        loc      int     location of vertex in array

    """
    loc = (N + coordinate[1]) * (2 * N + 1) + N + coordinate[0]
    return loc


def isLeaf_edge(N, G, edge):
    """
    isLeaf_edge checks if an edge is a leaf
    INPUT
        N        int        size parameter of box
        G        ig.Graph   graph in which to check for leaf
        edge     tuple      tuple with vertices (coord) of an edge

    OUTPUT
                 bool       True if edge is leaf in G
    """
    if(not G.are_connected(coord_to_array(N, edge[0]),
                           coord_to_array(N, edge[1]))):
        return print('ERROR edge ', edge, ' is not in graph')

    if(G.degree(coord_to_array(N, edge[0])) == 1
       or G.degree(coord_to_array(N, edge[1])) == 1):
        return True
    else:
        return False


def isLeaf_edge_id(G, edge):
    """
    isLeaf_edge_id checks if an edge is a leaf
    INPUT
        G        ig.Graph   graph in which to check for leaf
        edge     edge id    edge id (iGraph) of edge

    OUTPUT
                 bool       True if edge is leaf in G
    """
    if(G.degree(edge.source) == 1
       or G.degree(edge.target) == 1):
        return True
    else:
        return False
