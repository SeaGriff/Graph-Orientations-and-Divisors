from sage.graphs.graph import Graph
from sage.graphs.orientations import random_orientation

class CycleCocycleSystem(Graph):
    def __init__(self, inputGraph, _base_orientation = None, _base_edge = None, _base_vertex = None):
        assert type(inputGraph) == sage.graphs.graph.Graph, "Input is not a graph."
        assert not inputGraph.is_directed(), "Graph is directed."
        assert inputGraph.is_biconnected(), "Graph is not 2-edge connected."
        ad = inputGraph.adjacency_matrix()
        Graph.__init__(self, ad)
        self._pic = Sandpile(self)
        if _base_orientation == None:
            self._base_orientation = self.random_orientation()
        _big_theta_div = None
        _big_theta_orientation = None

    def set_big_theta_div(self, div):
        self._big_theta_div = div_op

    def set_big_theta_orientation(self, orientation):
        _big_theta_orientation = orientation

    def linking_add(self, A, B):
        assert self._big_theta_div != None, "Must set big theta divisor first."
        return A + B - self._big_theta_div

    def linking_orientation_add(self, U, W):
        return self.linear_orientation_class(self.linking_add(self.chern_class(A), self.chern_class(B)))

    # returns the underlying graph.
    def underlying_graph(self):
        return Graph(self.adjacency_matrix())

    def _repr_(self):
        return "A graph with a cycle-cocycle reversal system"

    # returns the Chern class of an orientation.
    def chern_class(self, orientation):
        D_add = SandpileDivisor(self._pic, {v:orientation.in_degree(v) for v in self.vertex_iterator()})
        return D_add - self._pic.all_k_div(1)

    # Takes O(D) of a divisor (currently requires deg D = g-1)
    def linear_orientation_class(self, div):
        assert div.deg() == self._pic.genus() - 1, "Currrently can only take O(D) for deg D = g-1."
        act_by = div - self.chern_class(self.base_orientation())
        return self.pic_0_action(self._base_orientation, act_by)

    # returns the base orientation
    def base_orientation(self):
        return self._base_orientation.copy()

    # performs the orientation equivalent of passing to the q-reduced representative of a divisor
    def q_red_orientation(self, orientation, q):
        self.make_paths(orientation, q)

    # flips oriented cuts until either every vertex is accessible by an oriented make_paths
    # from q, or if a target vertex is selected until the target is accessible.
    def make_paths(self, orientation, origin, target = None):
        vertex_set = set(self.vertices())
        reachable = reachable_vertices(orientation, origin)
        new_orientation = orientation.copy()
        while not target in reachable and not reachable == vertex_set:
            complement = vertex_set - reachable
            new_orientation.reverse_edges(new_orientation.edge_boundary(complement), multiedges=True)
            reachable = reachable_vertices(new_orientation, origin)
        return new_orientation

    # Applies a divisor of degree 0 to an orientation according to the canonical torsor action
    def pic_0_action(self, orientation, div):
        assert div.deg() == 0, "Only divisors of degree zero can act on orientations."
        D = SandpileDivisor(self._pic, div.copy())
        new_orientation = orientation.copy()
        D_pos = div_pos(self._pic, D, False)
        D_neg = div_pos(self._pic, -D, False)
        for i in range(len(D_pos)):
            new_orientation = self._gen_action(new_orientation, D_neg[i], D_pos[i])
        return new_orientation

    # Applies a divisor specifically of the form p - q to an orientation according
    # to the canonical torsor action
    def _gen_action(self, orientation, minusVertex, plusVertex):
        new_orientation = self.make_paths(orientation, plusVertex, minusVertex)
        P = next(new_orientation.shortest_simple_paths(plusVertex, minusVertex, report_edges = True))
        new_orientation.reverse_edges(P, multiedges = True)
        return new_orientation

    # Performs the divisor equivalent of flipping all edges in an orientation
    def div_op(self, div):
        return self._pic.canonical_divisor() - div

    # Returns a list of all theta characteristic divisors for G
    def theta_char_divisors(self):
        return [D for D in self._pic.picard_representatives(self._pic.genus() - 1) if D.is_linearly_equivalent(self.div_op(D))]

    # Returns a list of all theta characteristic divisors for G
    def theta_char_orientations(self):
        return [self.linear_orientation_class(D) for D in self.theta_char_divisors()]

    # Returns a list of representatives for the cycle cocycle system
    def orientation_representatives(self):
        return [self.linear_orientation_class(D) for D in self._pic.picard_representatives(self._pic.genus() - 1)]

    # iterates through deg g-1 picard_representatives and returns the first theta char found.
    # only intended to be used in the case there is a unique such.
    # if there is no theta char, returns False
    def big_theta_divisor(self):
        for D in self._pic.picard_representatives(self._pic.genus() - 1):
            if self.is_theta_char_divisor(D):
                return D
        return False

    def big_theta_orientation(self):
        result = self.big_theta_char_divisor()
        if result == False:
            return False
        return self.linear_orientation_class(result)

    def is_theta_char_divisor(self, div):
        return div.is_linearly_equivalent(self.div_op(div))

def div_pos(pic, div, dict_format = True):
    if dict_format:
        return SandpileDivisor(pic, {v: max(0, div[v]) for v in div.keys()})
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
