import random
import numpy as np
from scipy.stats import bernoulli
import ig_percoSimIterationAux as itaux

"""
this file contains wulffIteration_edges which compute one step of the Markov
chain
"""


def wulffIteration_edges(G, boundary_edges, interior_edges, p, N):
    """
    wulffIteration_edges perform one iteration step of the Markov chain
    given some state G

    INPUT
        G                   ig.Graph()  the state from which we compute the
                                        next step
        boundary_edges      list        list containing the boundary edges of
                                        the cluster from which we transition
        interior_edges      list        list containing the interior edges of
                                        the cluster from which we transition
        p                   float       number corresponding to the density for
                                        the percolation configuration
        N                   int         size of the box in which simulation
                                        is performed


    OUTPUT
        G                   ig.Graph()  updated graph
        boundary_edges      list        potentially updated boundary edges
        interior_edges      list        potentially updated interior edges
        N_C                 int         number of neighbourd of initial state C
        N_D                 int         number of neighbourd of next state D
        transitionProba     float       probability how likely it is to go from
                                        C to D (reversible)
        transitionProba_non float       probability how likely it is to go from
                                        C to D (not reversible, conjectured)
    """

    interior_edges = list(set(interior_edges))
    boundary_edges = list(set(boundary_edges))

    e_interior_frozen = random.choice(interior_edges)  # e_interior is frzn set
    e_interior = [v for v in e_interior_frozen]
    e_interior_ig = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                              itaux.coord_to_array(N, e_interior[1]))

    e_boundary_frozen = random.choice(boundary_edges)
    e_boundary = [v for v in e_boundary_frozen]

    edges_interior_updated = interior_edges[:]
    edges_boundary_updated = boundary_edges[:]

    edges_interior_updated.append(e_boundary_frozen)
    F = G.copy()
    # remove the boundary edge from the boundary list (is now in interior)
    for i in [0, 1]:
        try:
            edges_boundary_updated.remove(frozenset((e_boundary[i],
                                                     e_boundary[1 - i])))
        except ValueError:
            pass
    edges_boundary_updated.append(e_interior_frozen)

    G.add_edge(itaux.coord_to_array(N, e_boundary[0]),
               itaux.coord_to_array(N, e_boundary[1]))

    # remove the interior edge from the interioir list (is now in boundary)
    for i in [0, 1]:
        try:
            edges_interior_updated.remove(frozenset((e_interior[i],
                                                     e_interior[1 - i])))
        except ValueError:
            pass

    bridges = G.bridges()
    e_boundary_ig = G.get_eid(itaux.coord_to_array(N, e_boundary[0]),
                              itaux.coord_to_array(N, e_boundary[1]))

    # check the connectivite after adding the boundary edge
    if(e_interior_ig in bridges):
        if(not itaux.isLeaf_edge(N, G, e_interior)):
            G.delete_edges([e_boundary_ig])
            return G, boundary_edges, interior_edges, np.nan, np.nan, np.nan, np.nan

    e_interior_ig_updated = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                                      itaux.coord_to_array(N, e_interior[1]))
    incident_edges = G.incident(itaux.coord_to_array(N, (0, 0)))

    # check that the vertex at (0, 0) remains in the cluster
    if(e_interior_ig_updated in incident_edges
       and len(list(incident_edges)) <= 1):
        G.delete_edges([e_boundary_ig])  # must not fail silently
        return G, boundary_edges, interior_edges, np.nan, np.nan, np.nan, np.nan

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
                edges_boundary_updated.append(frozenset((itaux.array_to_coord(N, v),
                                                         itaux.array_to_coord(N, w))))

    # updating boundary by removing old boundary edges not incident to cluster
    # anymore
    e_interior_ig = G.get_eid(itaux.coord_to_array(N, e_interior[0]),
                              itaux.coord_to_array(N, e_interior[1]))
    e = G.es[e_interior_ig]
    V = [e.source, e.target]
    G.delete_edges([e_interior_ig])

    for v in V:
        if(len(list(G.incident(v))) == 0):
            nbrs_v = itaux.neighbours_vid(N, v)
            for w in nbrs_v:
                if(len(list(G.incident(w))) == 0):
                    try:
                        edges_boundary_updated.remove(frozenset((itaux.array_to_coord(N, v),
                                                                 itaux.array_to_coord(N, w))))
                    except ValueError:
                        pass
                    try:
                        edges_boundary_updated.remove(frozenset((itaux.array_to_coord(N, w),
                                                                 itaux.array_to_coord(N, v))))
                    except ValueError:
                        pass

    edges_removed.append((e_interior[0], e_interior[1]))

    delta_boundary = len(edges_boundary_updated) - len(boundary_edges)

    transition_rate = itaux.transition_rate_edges(p, delta_boundary)
    # transitionProba = transition_rate / (1 + transition_rate)
    transitionProba_non = min(transition_rate, 1)
    transitionProba, N_C, N_D = itaux.transition_proba_rev(F, p, N, boundary_edges,
                                                           edges_boundary_updated,
                                                           e_boundary, e_interior)

    if(bernoulli.rvs(transitionProba) == 1):
        edges_boundary_updated = list(set(edges_boundary_updated))
        edges_interior_updated = list(set(edges_interior_updated))
        """
        return G, edges_boundary_updated, edges_interior_updated, np.nan, np.nan, np.nan, np.nan
        """
        return G, edges_boundary_updated, edges_interior_updated, N_C, N_D, transitionProba, transitionProba_non

    else:
        G.add_edge(itaux.coord_to_array(N, e_interior[0]),
                   itaux.coord_to_array(N, e_interior[1]))
        e_bdry_eid = G.get_eid(itaux.coord_to_array(N, e_boundary[0]),
                               itaux.coord_to_array(N, e_boundary[1]))
        G.delete_edges([e_bdry_eid])
        """
        return G, boundary_edges, interior_edges, np.nan, np.nan, np.nan, np.nan
        """
        return G, boundary_edges, interior_edges, N_C, N_D, transitionProba, transitionProba_non