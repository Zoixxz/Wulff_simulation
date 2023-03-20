import random
from scipy.stats import bernoulli
import ig_percoSimIterationAux as itaux


def wulffIteration_edges(G, boundary_edges, interior_edges, p, N):

    interior_edges = list(set(interior_edges))
    boundary_edges = list(set(boundary_edges))

    e_interior = random.choice(interior_edges)
    e_interior_ig = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                              itaux.coord_to_array(N, e_interior[1]))

    e_boundary = random.choice(boundary_edges)

    edges_interior_updated = interior_edges[:]
    edges_boundary_updated = boundary_edges[:]

    edges_interior_updated.append(e_boundary)
    F = G.copy()
    # remove the boundary edge from the boundary list (is now in interior)
    for i in [0, 1]:
        try:
            edges_boundary_updated.remove((e_boundary[i], e_boundary[1 - i]))
        except ValueError:
            pass
    edges_boundary_updated.append(e_interior)

    G.add_edge(itaux.coord_to_array(N, e_boundary[0]),
               itaux.coord_to_array(N, e_boundary[1]))

    # remove the interior edge from the interioir list (is now in boundary)
    for i in [0, 1]:
        try:
            edges_interior_updated.remove((e_interior[i], e_interior[1 - i]))
        except ValueError:
            pass

    bridges = G.bridges()
    e_boundary_ig = G.get_eid(itaux.coord_to_array(N, e_boundary[0]),
                              itaux.coord_to_array(N, e_boundary[1]))

    # check the connectivite after adding the boundary edge
    if(e_interior_ig in bridges):
        if(not itaux.isLeaf_edge(N, G, e_interior)):
            G.delete_edges([e_boundary_ig])
            return G, boundary_edges, interior_edges

    e_interior_ig_updated = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                                      itaux.coord_to_array(N, e_interior[1]))
    incident_edges = G.incident(itaux.coord_to_array(N, (0, 0)))

    # check that the vertex at (0, 0) remains in the cluster
    if(e_interior_ig_updated in incident_edges
       and len(list(incident_edges)) <= 1):
        G.delete_edges([e_boundary_ig])  # must not fail silently
        return G, boundary_edges, interior_edges

    # update the boundary and graph
    edges_removed = []
    # updating boundary by adding possibly new boundary edges
    e_boundary_ig = G.get_eid(itaux.coord_to_array(N, e_boundary[0]),
                              itaux.coord_to_array(N, e_boundary[1]))
    e = G.es[e_boundary_ig]
    V = [e.source, e.target]
    for v in V:
        nbrs_v = itaux.neighbours_vid(N, v)
        for w in nbrs_v:
            if(not G.are_connected(v, w)):
                edges_boundary_updated.append((itaux.array_to_coord(N, v),
                                               itaux.array_to_coord(N, w)))

    # updating boundary by removing old boundary edges not indicent to cluster
    # anymore
    e_interior_ig = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                              itaux.coord_to_array(N, e_interior[1]))
    e = G.es[e_interior_ig]
    V = [e.source, e.target]
    for v in V:
        nbrs_v = itaux.neighbours_vid(N, v)
        for w in nbrs_v:
            if(not G.are_connected(v, w)):
                try:
                    edges_boundary_updated.remove((itaux.array_to_coord(N, v),
                                                   itaux.array_to_coord(N, w)))
                except ValueError:
                    pass
                try:
                    edges_boundary_updated.remove((itaux.array_to_coord(N, w),
                                                   itaux.array_to_coord(N, v)))
                except ValueError:
                    pass

    edges_removed.append((e_interior[0], e_interior[1]))
    G.delete_edges([e_interior_ig])

    delta_boundary = -(len(edges_boundary_updated) - len(boundary_edges))

    transition_rate = itaux.transition_rate_edges(p, delta_boundary)
    # transitionProba = transition_rate / (1 + transition_rate)
    transitionProba = min(transition_rate, 1)
    """
    transitionProba_2 = itaux.transition_proba_rev(F, p, N, interior_edges,
                                                   boundary_edges,
                                                   edges_interior_updated,
                                                   edges_boundary_updated,
                                                   e_boundary, e_interior)
    """
    if(bernoulli.rvs(transitionProba) == 1):
        return G, edges_boundary_updated, edges_interior_updated
    else:
        G.add_edge(itaux.coord_to_array(N, e_interior[0]),
                   itaux.coord_to_array(N, e_interior[1]))
        e_bdry_eid = G.get_eid(itaux.coord_to_array(N, e_boundary[0]),
                               itaux.coord_to_array(N, e_boundary[1]))
        G.delete_edges([e_bdry_eid])
        return G, boundary_edges, interior_edges
