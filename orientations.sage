from sage.graphs.graph import Graph
from sage.graphs.orientations import random_orientation

class CycleCocycleSystem(Graph):
    def __init__(self, inputGraph, _pic = None, _base_orientation = None):
        assert type(inputGraph) == sage.graphs.graph.Graph, "Input is not a graph."
        assert not inputGraph.is_directed(), "Graph is directed."
        assert inputGraph.is_biconnected(), "Graph is not 2-edge connected."
        ad = inputGraph.adjacency_matrix()
        Graph.__init__(self, ad)
        if _pic == None:
            self._pic = Sandpile(self)
        if _base_orientation == None:
            self._base_orientation = self.random_orientation()

    def picard_group(self):
        return self._pic.copy()

    def _repr_(self):
        return "A graph with a cycle-cocycle reversal system."

    def chern_class(self, orientation):
        D_add = SandpileDivisor(self._pic, {v:orientation.in_degree(v) for v in self.vertex_iterator()})
        return D_add - self._pic.all_k_div(1)

    def base_orientation(self):
        return self._base_orientation.copy()

    def q_red(self, orientation, q):
        vertex_set = set(self.vertices())
        reachable = reachable_vertices(orientation, q)
        new_orientation = orientation.copy()
        while not reachable == vertex_set:
            complement = vertex_set - reachable
            new_orientation.reverse_edges(new_orientation.edge_boundary(complement), multiedges=True)
            reachable = reachable_vertices(new_orientation, q)
        return new_orientation

    def gen_action(self, orientation, minusVertex, plusVertex):
        new_orientation = self.q_red(orientation, plusVertex)
        P = next(new_orientation.shortest_simple_paths(plusVertex, minusVertex, report_edges = True))
        new_orientation.reverse_edges(P, multiedges = True)
        return new_orientation

def reachable_vertices(G, q):
    reachable = {q}
    reachable_it(G, q, reachable)
    return reachable

def reachable_it(G, v, reachable):
    reachable.add(v)
    for w in G.vertex_boundary(reachable):
        reachable_it(G, w, reachable)
