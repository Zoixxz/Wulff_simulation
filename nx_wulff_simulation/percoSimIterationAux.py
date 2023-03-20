def neighbours(x):
    """
    takes x [tuple] a vertex as input and return
    a list of all neighbours in Z^2
    """
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
            for i in range(len(iterate_neighbours))]


def transition_rate_edges(p, Delta_boundary):
    return (1 - p)**Delta_boundary


def isLeaf(G, edge):
    # check if the edge is in the graph
    if(not G.has_edge(edge[0], edge[1])):
        return print('ERROR edge ', edge, ' is not in graph')

    if(G.degree(edge[0]) == 1 or G.degree(edge[1]) == 1):
        return True
    else:
        return False
