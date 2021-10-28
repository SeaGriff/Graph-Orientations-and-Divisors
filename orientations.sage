from sage.graphs.graph import Graph
from sage.graphs.orientations import random_orientation

class CycleCocycleSystem(Graph):
    def __init__(self, inputGraph, _pic = None, _base_orientation = None, _base_edge = None, _base_vertex = None):
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
        return self._pic

    def _repr_(self):
        return "A graph with a cycle-cocycle reversal system"

    def chern_class(self, orientation):
        D_add = SandpileDivisor(self._pic, {v:orientation.in_degree(v) for v in self.vertex_iterator()})
        return D_add - self._pic.all_k_div(1)

    def linear_orientation_class(self, div):
        assert div.deg() == self._pic.genus() - 1, "Currrently can only take O(D) for deg D = g-1."
        act_by = div - self.chern_class(self.base_orientation())
        return self.pic_0_action(self._base_orientation, act_by)

    def base_orientation(self):
        return self._base_orientation.copy()

    def q_red_orientation(self, orientation, q):
        self.make_paths(orientation, q)

    def make_paths(self, orientation, origin, target = None):
        vertex_set = set(self.vertices())
        reachable = reachable_vertices(orientation, origin)
        new_orientation = orientation.copy()
        while not target in reachable and not reachable == vertex_set:
            complement = vertex_set - reachable
            new_orientation.reverse_edges(new_orientation.edge_boundary(complement), multiedges=True)
            reachable = reachable_vertices(new_orientation, origin)
        return new_orientation

    def pic_0_action(self, orientation, div):
        assert div.deg() == 0, "Only divisors of degree zero can act on orientations."
        D = SandpileDivisor(self._pic, div.copy())
        new_orientation = orientation.copy()
        D_pos = get_pos_part(D, False)
        D_neg = get_pos_part(-D, False)
        for i in range(len(D_pos)):
            new_orientation = self._gen_action(new_orientation, D_neg[i], D_pos[i])
        return new_orientation

    def _gen_action(self, orientation, minusVertex, plusVertex):
        new_orientation = self.make_paths(orientation, plusVertex, minusVertex)
        P = next(new_orientation.shortest_simple_paths(plusVertex, minusVertex, report_edges = True))
        new_orientation.reverse_edges(P, multiedges = True)
        return new_orientation

    def div_op(self, div):
        return self._pic.canonical_divisor() - div

    def get_theta_char_divisors(self):
        return [D for D in self._pic.picard_representatives(self._pic.genus() - 1) if D.is_linearly_equivalent(self.div_op(D))]

    def get_theta_char_orientations(self):
        return [self.linear_orientation_class(D) for D in self.get_theta_char_divisors()]

def get_pos_part(div, dict_format = True):
    if dict_format:
        return SandpileDivisor(self._pic, {v: max(0, div[v]) for v in div.keys()})
    output_list = []
    for v in div.keys():
        if div[v] > 0:
            output_list.extend([v] * int(div[v]))
    return output_list

def reachable_vertices(G, q):
    reachable = {q}
    reachable_it(G, q, reachable)
    return reachable

def reachable_it(G, v, reachable):
    reachable.add(v)
    for w in G.vertex_boundary(reachable):
        reachable_it(G, w, reachable)
