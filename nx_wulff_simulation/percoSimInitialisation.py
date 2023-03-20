from scipy.stats import bernoulli


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


def neighbours(x):
    """
    takes x [tuple] a vertex as input and 
    returns a list of all neighbours in Z^2
    """
    iterate_neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return [tuple(map(lambda i, j: i + j, x, iterate_neighbours[i]))
            for i in range(len(iterate_neighbours))]


def graphDistanceMax(v):
    # v: vertex with coordinates
    # returns max norm of the vertes/coordinate
    # works is any dimension
    dist = max([abs(v[i]) for i in range(len(v))])
    return dist


def graphDistanceOne(v):
    """
    computes the one norm of a vector v
    """
    dist = sum([abs(v[i]) for i in range(len(v))])
    return dist


def crateInitialConfiguration(G, N, m, p):
    """
    function takes as inputs an nx graph G and a size 
    for the Graph to be created
    returns the initial configuration of the graph which we wan
    Create the correct vertices of the graph from -N to N in a square grid

    intputs:
    G	intitialized nx.Graph(), must be empty
    N 	size of vertex box for simulation
    m 	size of box of initial configuration
    p 	probabilistic parameter


    returns:
    G 					[an nx graph] 	a graph with the correct connectivity
    V 					[array] 		of vertex coordinates of all vertices in the graph G
    exteriorBoundary 	[array]			of vertex coordinates exterior boundary of C_0
    cluster0			[array]			of vertex coordinates of the cluster containing 0

    """

    # Create the graph

    V = []  # empty array of Vertices
    # build the vertex set [a list of vertices]
    # V is a 2N+1 time 2N+1 array
    # V[0] = [-N, -N] and V[1] = [-N+1, -N]

    for j in range(-N, N + 1):
        for i in range(-N, N + 1):
            V.append([i, j])
    V = [tuple(V[i]) for i in range(len(V))]

    # add all the vertices/nodes to G
    for i in range(len(V)):
        # need to convert to tuple as list is not iterable object
        G.add_node(tuple(V[i]), pos=tuple(V[i]))

    # add all the eges (random) to G
    # iterate over horizontal (x-dir) edges

    for i in range((2 * N + 1) * (2 * N + 1) - 1):
        if(isNeighbour(V[i], V[i + 1]) and bernoulli.rvs(p) == 1):
            G.add_edge(V[i], V[i + 1])

    # iterate over vertical (y-dir) edges
    for i in range((2 * N + 1) * 2 * N):
        if(isNeighbour(V[i], V[i + (2 * N + 1)]) and bernoulli.rvs(p) == 1):
            G.add_edge(V[i], V[i + (2 * N + 1)])

    for v in V:
        if(graphDistanceOne(v) <= m):
            neighbours_v = neighbours(v)
            for w in neighbours_v:
                if(graphDistanceOne(v) > m):
                    if(G.has_edge(v, m)):
                        G.remove_edge(v, m)
                else:
                    G.add_edge(v, m)

    exteriorBoundary = [v for v in V if (graphDistanceOne(
        v) == m + 1 and min([graphDistanceOne(w) for
                             w in neighbours(v)]) <= m)]
    cluster0 = [v for v in V if graphDistanceOne(v) <= m]

    return G, V, exteriorBoundary, cluster0
