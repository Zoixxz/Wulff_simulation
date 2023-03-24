import random
import igraph as ig
from scipy.stats import bernoulli
import ig_percoSimIterationAux as itaux


def wulffIteration_edges(G, boundary_vertices, interior_vertices, p, N):

    interior_vertices = list(set(interior_vertices))
    boundary_vertices = list(set(boundary_vertices))

    v_interior = random.choice(interior_vertices)
    v_interior_ig = itaux.coord_to_array(N, v_interior)

    v_boundary = random.choice(boundary_vertices)

    vertices_interior_updated = interior_vertices[:]
    vertices_boundary_updated = boundary_vertices[:]

    vertices_interior_updated.append(v_boundary)
    # remove the boundary edge from the boundary list (is now in interior)
    try:
        vertices_boundary_updated.remove(v_boundary)
    except ValueError:
        pass

    vertices_boundary_updated.append(v_interior)

    # add the edges to the graph which are in the cluster and connecting
    # the boundary vertex (as is approprate for siter perco)
    neighbours_v_boundary = itaux.neighbours(v_boundary)
    edges_added = []
    # remove the interior edge from the interior list (is now in boundary)
    for w in neighbours_v_boundary:
        if w in interior_vertices:
            G.add_edge(itaux.coord_to_array(N, w),
                       itaux.coord_to_array(N, v_boundary))
            edges_added.append((itaux.coord_to_array(N, w),
                                itaux.coord_to_array(N, v_boundary)))

    vertices_interior_updated.remove(v_interior)

    articulation_points = G.articulation_points()

    # check the connectivite after adding the boundary edge
    if(v_interior_ig in articulation_points):
        G.delete_edges(edges_added)
        return G, boundary_vertices, interior_vertices

    e_interior_ig_updated = G.incident(itaux.coord_to_array(N, v_interior))
    incident_edges = G.incident(itaux.coord_to_array(N, (0, 0)))

    # check that the vertex at (0, 0) remains in the cluster
    if(any(e in incident_edges) for e in e_interior_ig_updated
       and len(list(incident_edges)) <= 1):
        G.delete_edges(edges_added)  # must not fail silently
        return G, boundary_vertices, interior_vertices

    # update the boundary and graph
    edges_removed = []
    # updating boundary by adding possibly new boundary edges
    v_boundary_ig = G.get_vid(itaux.coord_to_array(N, v_boundary))
    v = G.vs[v_boundary_ig]
    V = [e.source, e.target]
    for v in V:
        nbrs_v = itaux.neighbours_vid(N, v)
        for w in nbrs_v:
            if(not G.are_connected(v, w)):
                edges_boundary_updated.append((itaux.array_to_coord(N, v),
                                               itaux.array_to_coord(N, w)))
    for w in itaux.neigbours(v_boundary):
        if(not w in ig.neighbours(G, v_boundary_ig)):
            vertices_boundary_updated.append(itaux.array_to_coord(N, w))

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
