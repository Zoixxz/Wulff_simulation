import random
from scipy.stats import bernoulli
import networkx as nx
import percoSimIterationAux as itaux
import time as t


def neighbours(x):
    """
    takes x [tuple] a vertex as input and return a list
    of all neighbours in Z^2
    """
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
            for i in range(len(iterate_neighbours))]


def wulffIteration_edges(G, boundary_edges, interior_edges, p):

    interior_edges = list(set(interior_edges))
    boundary_edges = list(set(boundary_edges))

    e_interior = random.choice(interior_edges)
    e_boundary = random.choice(boundary_edges)

    edges_interior_updated = interior_edges[:]
    edges_boundary_updated = boundary_edges[:]

    edges_interior_updated.append(e_boundary)

    for i in [0, 1]:
        try:
            edges_boundary_updated.remove(tuple([e_boundary[i],
                                                 e_boundary[1 - i]]))
        except ValueError:
            pass
    edges_boundary_updated.append(e_interior)

    G.add_edge(e_boundary[0], e_boundary[1])

    Cluster_new_connectivity = G.edge_subgraph(edges_interior_updated)

    Cluster_old = G.edge_subgraph(interior_edges)

    for i in [0, 1]:
        try:
            edges_interior_updated.remove(tuple([e_interior[i],
                                                 e_interior[1 - i]]))
        except ValueError:
            pass

    Cluster_new = G.edge_subgraph(edges_interior_updated)

    # check in e_interior is the only edge connecting to 0
    bridges = list(nx.bridges(Cluster_new_connectivity))

    if not Cluster_new.has_node((0, 0)):
        G.remove_edge(e_boundary[0], e_boundary[1])  # must not fail silently
        return G, boundary_edges, interior_edges, list(bridges)

    # check connectivity
    if(e_interior in bridges or tuple(reversed(e_interior)) in bridges):
        if(not itaux.isLeaf(Cluster_new_connectivity, e_interior)):
            G.remove_edge(e_boundary[0], e_boundary[1])  # not fail silently
            return G, boundary_edges, interior_edges, list(bridges)

    # update the boundary and graph
    edges_removed = []

    for i in [0, 1]:
        if(not Cluster_old.has_node(e_boundary[i])):
            nbrs_e_0 = neighbours(e_boundary[i])
            for v in nbrs_e_0:
                if(not Cluster_new.has_edge(v, e_boundary[i])):
                    edges_boundary_updated.append((v, e_boundary[i]))
                    # possibly have added duplicates
                    if(G.has_edge(v, e_boundary[i])):
                        edges_removed.append((v, e_boundary[i]))
                        G.remove_edge(v, e_boundary[i])
    for i in [0, 1]:
        if not Cluster_new.has_node(e_interior[i]):
            nbrs_e_0 = neighbours(e_interior[i])
            for v in nbrs_e_0:
                if(not Cluster_new.has_node(v)):
                    try:
                        edges_boundary_updated.remove(tuple([v,
                                                             e_interior[i]]))
                    except ValueError:
                        pass
                    try:
                        edges_boundary_updated.remove(tuple([e_interior[i],
                                                             v]))
                    except ValueError:
                        pass

    if(G.has_edge(e_interior[0], e_interior[1])):
        edges_removed.append((e_interior[0], e_interior[1]))
        G.remove_edge(e_interior[0], e_interior[1])

    delta_boundary = len(edges_boundary_updated) - len(boundary_edges)

    transition_rate = itaux.transition_rate_edges(p, delta_boundary)
    transitionProba = transition_rate / (1 + transition_rate)
    # print('transition_rate ', transition_rate)
    # print('transitionProba ', transitionProba)
    if(bernoulli.rvs(transitionProba) == 1):
        return G, edges_boundary_updated, edges_interior_updated, list(bridges)
    else:
        # we reverte the changes in the graph
        G.remove_edge(e_boundary[0], e_boundary[1])  # must not fail silently
        for e in edges_removed:
            G.add_edge(e[0], e[1])
        return G, boundary_edges, interior_edges, list(bridges)
