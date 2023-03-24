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
    takes x [tuple] a vertex as input and return
    a list of all neighbours in Z^2 with the vertex coordinates
    """
    x = array_to_coord(N, y)
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    V = [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
         for i in range(len(iterate_neighbours))]
    return [coord_to_array(N, v) for v in V]


def transition_rate_edges(p, Delta_boundary):
    return (1 - p) ** Delta_boundary


def transition_proba_rev(G, p, N, interior_edges, boundary_edges,
                         interior_edges_updated, boundary_edges_updated,
                         e_bdry, e_interior):
    """
    given:
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
        TP                      transition probability to the other state
    """
    edge_count = G.ecount()
    F = G.copy()
    N_C = 0
    N_D = 0
    # compute N(C), the neighbouring states of C for the relation R
    for v in boundary_edges:
        F.add_edge(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        bridges = F.es[F.bridges()]
        bridges_number = len([e for e in bridges
                              if not isLeaf_edge_id(F, e)])
        N_C += edge_count - bridges_number
        e_bdry_id = F.get_eid(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        F.delete_edges([e_bdry_id])

    # compute N(D), the neighbouring states of D for the relation R
    # update the graph F so that we have configuaration D
    F.add_edge(coord_to_array(N, e_bdry[0]), coord_to_array(N, e_bdry[1]))
    e_id = F.get_eid(coord_to_array(N, e_interior[0]),
                     coord_to_array(N, e_interior[1]))
    F.delete_edges([e_id])

    for v in boundary_edges_updated:
        F.add_edge(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        bridges = F.es(F.bridges())
        bridges_number = len([e for e in bridges
                              if not isLeaf_edge_id(F, e)])
        N_D += edge_count - bridges_number
        e_bdry_id = F.get_eid(coord_to_array(N, v[0]), coord_to_array(N, v[1]))
        F.delete_edges([e_bdry_id])
    Delta_boundary = len(boundary_edges_updated) - len(boundary_edges)
    return min((1 - p) ** Delta_boundary * (N_C / N_D), 1)


def array_to_coord(N, location):
    # igraph counts vertices starting at 0 (in python)
    m_y, m_x = divmod(location, 2 * N + 1)
    return (m_x - N, m_y - N)


def coord_to_array(N, coordinate):
    # igraph counts vertices starting at 0 (in python)
    loc = (N + coordinate[1]) * (2 * N + 1) + N + coordinate[0]
    return loc


def isLeaf_edge(N, G, edge):
    # check if the edge is in the graph
    if(not G.are_connected(coord_to_array(N, edge[0]),
                           coord_to_array(N, edge[1]))):
        return print('ERROR edge ', edge, ' is not in graph')

    if(G.degree(coord_to_array(N, edge[0])) == 1
       or G.degree(coord_to_array(N, edge[1])) == 1):
        return True
    else:
        return False


def isLeaf_edge_id(G, edge):
    if(G.degree(edge.source) == 1
       or G.degree(edge.target) == 1):
        return True
    else:
        return False
