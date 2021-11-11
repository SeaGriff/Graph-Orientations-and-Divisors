import itertools
from sage.graphs.generic_graph import GenericGraph
from sage.sandpiles.sandpile import SandpileDivisor
from sage.graphs.digraph import DiGraph

# SandpileDivisor Methods

def div_op(div):
    """ Performs the divisor equivalent of flipping all
    edges in an orientation (including biorienting <-> unorienting) """
    return div.sandpile().canonical_divisor() - div


def div_pos(div, divisor_format=True):
    """ Returns the positive part of a divisor, or if divisor_format is False,
    a list of all positive entries (with multiplicity) in the
    divisor. """
    if divisor_format:
        return SandpileDivisor(div.sandpile(),
                               {v: max(0, div[v]) for v in div.keys()})
    output_list = []
    for v in div.keys():
        if div[v] > 0:
            output_list.extend([v] * int(div[v]))
    return output_list


# DiGraph methods

def reachable_from_vertex(G, q):
    """ Returns the set of vertices reachable by oriented paths from q,
    a vertex in a DiGraph. """
    return reachable_from_vertices({q})

def reachable_from_vertices(diG, X):
    """ Returns the set of vertices reachable by oriented paths from X,
    a collection of vertices in a DiGraph. """
    V = set(X)
    to_add = diG.edge_boundary(V)
    while len(to_add) != 0:
        X.union({e[1] for e in to_add})
        to_add = diG.edge_boundary(V)
    return V

# Generic graph methods

def induces_connected_subgraph(G, vertices):
    """ Returns whether a set of vertices determines a connected induced
    subgraph. Treats the empty graph as not connected. """
    if len(vertices) == 0:
        return False
    return G.subgraph(vertices).is_connected()


def eulerian_bipartition(G):
    """ returns a set of vertices V of a graph G,
    # such that G[V] and G[V^c] are eulerian.
    # We have V = V(G) iff the graph has purely even degrees. """
    if G.has_multiple_edges():
        preprocess_G = Graph([G.vertices(), []])
        for (v, w) in itertools.combinations(G.vertices(), 2):
            if is_odd(len(G.edge_boundary([v], [w], labels=False))):
                preprocess_G.add_edge(v, w)
        result = eulerian_bipartition_recur(preprocess_G, ([], []))
    else:
        result = eulerian_bipartition_recur(G, ([], []))
    return result[0]


def eulerian_bipartition_recur(G, partition):
    """ The recursion for the eulerian bipartition algorithm. """
    for v in G.vertices():
        if is_odd(G.degree(v)):
            smallG = G.copy()
            smallG.delete_vertex(v)
            neighbs = G.neighbors(v)
            for vertex_pair in itertools.combinations(neighbs, 2):
                if G.has_edge(vertex_pair[0], vertex_pair[1]):
                    smallG.delete_edge(vertex_pair[0], vertex_pair[1])
                else:
                    smallG.add_edge(vertex_pair[0], vertex_pair[1])
            partition = eulerian_bipartition_recur(smallG, partition)
            if is_odd(len(set(partition[0]).intersection(G.neighbors(v)))):
                partition[1].extend([v])
            else:
                partition[0].extend([v])
            return partition
    return (G.vertices(), [])


def partition_to_theta_char_orientation(G, V, show=False):
    """ Accepts a graph with a nonempty subset of the vertices. The subset and
    # complement are assumed to induce two Eulerian subgraphs.
    # Returns an orientation which is Eulerian on the two subgraphs and
    # has a consistently oriented cut between them. """
    if set(V) == set(G.vertices()):
        return G.eulerian_orientation()
    V_comp = vertex_complement(G, V)
    G1 = G.subgraph(V).eulerian_orientation()
    G2 = G.subgraph(V_comp).eulerian_orientation()
    result = DiGraph([G.vertices(), []], multiedges=G.allows_multiple_edges())
    result.add_edges(G.edge_boundary(V))
    result.reverse_edges(result.incoming_edge_iterator(V))
    result.add_edges(G1.edges())
    result.add_edges(G2.edges())
    if show:
        result.show(vertex_colors={'b': V, 'r': V_comp})
    return result


def vertex_complement(G, V):
    """ Returns the complement of a subset of vertices. """
    return set(G.vertices()) - set(V)

SandpileDivisor.div_op = div_op
SandpileDivisor.div_pos = div_pos
DiGraph.reachable_from_vertex = reachable_from_vertex
DiGraph.reachable_from_vertices = reachable_from_vertices
GenericGraph.eulerian_bipartition = eulerian_bipartition
GenericGraph.vertex_complement = vertex_complement
GenericGraph._parti_to_or = partition_to_theta_char_orientation
