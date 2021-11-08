""" Builds a cycle-cocycle system on top of a graph and implements associated
methods, in particular those that relate orientations to divisors. """

import itertools
from sage.graphs.graph import Graph
from copy import copy

class SuperDiGraph(DiGraph):
    """ Accepts a DiGraph and implements methods to view it as a super
    directed graph (a graph where edges may additionally be unoriented, or
    unoriented in both directions) """
    def __init__(self, G, data, bi=set(), unori=set()):
        assert bi.isdisjoint(unori), "Edge cannot be both bioriented and unoriented."

        DiGraph.__init__(self, data)

        self._bi, self._unori = bi, unori
        if isinstance(G, CycleCocycleSystem):
            self._CCS = G
        else:
            self._CCS = CycleCocycleSystem(G)

    def __repr__(self):
        return ("A superorientation on a graph with {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def copy(self):
        return SuperDiGraph(self._CCS, self, self._bi, self._unori)

    def biori(self):
        return self._bi.copy()

    def unori(self):
        return self._unori.copy()

    def CCS(self):
        return self._CCS.copy()

    def set_unori(self, X):
        assert X.isdisjoint(self._bi), "Edge cannot be both bioriented and unoriented."
        self._unori = X

    def set_biori(self, X):
            assert X.isdisjoint(self._unori), "Edge cannot be both bioriented and unoriented."
            self._bi = X

    def op(self):
        return SuperDiGraph(self._CCS, self, self._unori, self._bi)

    def is_equivalent(self, U):
        return self.chern_class().is_linearly_equivalent(self.chern_class(W))

    def sources(self):
        return self._make_algo_graph().sources()

    def unfurl(self, X, check=True):
        U = self.copy()
        if check:
            assert len(self._bi) == 0, "Currently only supports pure partial orientations."
            assert len(X) != 0, "Input must be nonempty."
            assert induces_connected_subgraph(Graph(self), X), "Input must induce a connected subgraph."
        return self._unfurl_it(X, U)

    def _unfurl_it(self, X, U):
        if set(X).issubset(set(U.sources()):
            unori_U = Graph(U)
            X_boundary = unori_U.vertex_boundary(X)
            for v in X_boundary:
                
            return
        return U

    def _make_algo_graph(self):
        result = DiGraph(self)
        result.delete_edges(self._unori)
        result.add_edges({(e[1], e[0]) for e in self._bi})
        return result

    def chern_class(self):
        """ returns the Chern class of the orientation. """
        D = self._CCS.pic().all_k_div(-1)
        for e in self.edges():
            if e not in self._unori:
                D[e[1]] += 1
                if e in self._bi:
                    D[e[0]] += 1
        return D

    def is_theta_char(self):
        return self._CCS.is_theta_char(self)

class CycleCocycleSystem(Graph):
    """ Accepts a graph and returns a cycle-cocycle system, a graph
    with methods outputting various data related to the system. """
    def __init__(self, data, base_orientation=None, base_edge=None):
        # Build the graph
        Graph.__init__(self, data)

        # Error checking
        assert self.is_biconnected(), "Graph is not 2-edge connected."

        # Initialize internal variables
        self._big_theta_div = None
        self._big_theta_orientation = None

        # Build associated objects
        self._pic = Sandpile(self)
        if base_orientation is None:
            self._base_orientation = SuperDiGraph(self,
                                                  self.random_orientation())
        else:
            self._base_orientation = SuperDiGraph(self,
                                                  base_orientation)
        if base_edge is None:
            self._base_edge = self.edges()[0]

    def _repr_(self):
        return ("A graph with a cycle-cocycle reversal system, on {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def copy(self):
        return CycleCocycleSystem(self, self._base_orientation, self._base_edge)

    # Set or retrieve attached objects and internal variables

    def set_big_theta(self, T, check=True):
        """ Set the big theta divisor and orientation.
        Accepts a divisor or orientation. """
        if isinstance(T, SuperDiGraph):
            if check:
                assert self.is_theta_char(T), "Input must be a \
                theta characteristic orientation"
            self._big_theta_div = T.chern_class()
            self._big_theta_orientation = T
        else:
            if check:
                assert self.is_theta_char(T), "Input must be a theta \
                characteristic divisor"
            self._big_theta_div = T
            self._big_theta_orientation = self.linear_orientation_class(D)

    def set_base_edge(self, e):
        """ Sets the base edge """
        self._base_edge = e

    def set_base_orientation(self, G):
        """ Sets the base orientation """
        self._base_orientation = G

    def pic(self):
        """ return the picard group of the graph """
        return self._pic

    def base_orientation(self):
        """ returns the base orientation """
        return self._base_orientation.copy()

    def base_edge(self):
        """ returns the base edge """
        return self._base_edge

    def big_theta_orientation(self):
        """ returns the currently set big theta orientation """
        return self._big_theta_orientation.copy()

    def big_theta_divisor(self):
        """ returns the currently set big theta divisor """
        return self._big_theta_div

    def theta_char_divisors(self):
        """ Returns a list of all theta characteristic divisors for G. Slow """
        return [D for D in self._pic.picard_representatives(self._pic.genus() -
                1) if self.is_theta_char(D)]

    def theta_char_orientations(self):
        """ Returns a list of all theta characteristic orientations for G.
        Slow """
        return [self.linear_orientation_class(D) for D in
                self.theta_char_divisors()]

    def sample_theta_char_div(self):
        """ returns a single theta characteristic divisor. Fast """
        return self.sample_theta_char_orientation().chern_class()

    def sample_theta_char_orientation(self, show=False):
        """ returns a single theta characteristic orientation. Fast """
        G = partition_to_theta_char_orientation(Graph(self),
                        eulerian_bipartition(Graph(self), show))
        return SuperDiGraph(self, G)

    def orientation_representatives(self):
        """ Returns a list of representatives for the cycle cocycle system """
        return [self.linear_orientation_class(D) for D
                in self._pic.picard_representatives(self._pic.genus() - 1)]

    # Operations on orientations

    def q_red_orientation(self, orientation, q):
        """ performs the orientation equivalent of passing to the q-reduced
        representative of a divisor """
        self.make_paths(orientation, q)

    def make_paths(self, orientation, origin, target=None):
        """ flips oriented cuts til either every vertex is accessible by an
        oriented path from q, or, if a target vertex is selected, until
        the target is accessible. """
        vertex_set = set(self.vertices())
        reachable = reachable_vertices(orientation, origin)
        new_orientation = orientation.copy()
        while (target not in reachable) and not reachable == vertex_set:
            complement = vertex_set - reachable
            new_orientation.reverse_edges(
                new_orientation.edge_boundary(complement), multiedges=True)
            reachable = reachable_vertices(new_orientation, origin)
        return new_orientation

    def pic_0_action(self, orientation, div):
        """ Applies a divisor of degree 0 to an orientation according to the
        canonical torsor action """
        assert div.deg() == 0, "Only divisors of degree zero can act \
            on orientations."
        D = SandpileDivisor(self._pic, div.copy())
        new_ori = orientation.copy()
        D_pos = div_pos(self._pic, D, False)
        D_neg = div_pos(self._pic, -D, False)
        for i in range(len(D_pos)):
            new_ori = self._gen_action(new_ori, D_neg[i], D_pos[i])
        return new_ori

    def _gen_action(self, orientation, minusVertex, plusVertex):
        """ Applies a divisor specifically of the form p - q to an orientation
        according to the canonical torsor action """
        new_ori = self.make_paths(orientation, plusVertex, minusVertex)
        P = next(new_ori.shortest_simple_paths(plusVertex,
                                            minusVertex, report_edges=True))
        new_ori.reverse_edges(P, multiedges=True)
        return new_ori

    # Operations on divisors

    def linear_orientation_class(self, div):
        """ takes O(D) of a divisor (currently requires deg D = g-1) """
        assert div.deg() == self._pic.genus() - 1, "Currrently can only take \
        O(D) for deg D = g-1."
        act_by = div - self.base_orientation().chern_class()
        return SuperDiGraph(self,
                            self.pic_0_action(self._base_orientation, act_by))

    # Misc

    def is_theta_char(self, T):
        """  """
        if isinstance(T, SuperDiGraph):
            if len(T.biori_set()) != len(T.unori_set()):
                return False
            T = T.chern_class()
        return T.is_linearly_equivalent(div_op(self._pic, T))


def div_op(S, div):
    """ Performs the divisor equivalent of flipping all
    edges in an orientation """
    return S.canonical_divisor() - div


def div_pos(pic, div, dict_format=True):
    """  """
    if dict_format:
        return SandpileDivisor(pic, {v: max(0, div[v]) for v in div.keys()})
    output_list = []
    for v in div.keys():
        if div[v] > 0:
            output_list.extend([v] * int(div[v]))
    return output_list


def reachable_vertices(G, q):
    """  """
    reachable = {q}
    reachable_it(G, q, reachable)
    return reachable


def reachable_it(G, v, reachable):
    """  """
    reachable.add(v)
    for w in G.vertex_boundary(reachable):
        reachable_it(G, w, reachable)


def induces_connected_subgraph(G, vertices):
    if len(vertices) == 0:
        return False
    return G.subgraph(vertices).is_connected()


def eulerian_bipartition(G, show=False):
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
    if show:
        G.show(vertex_colors={'b': result[0], 'r': result[1]})
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
    V_complement = set(G.vertices()) - set(V)
    G1 = G.subgraph(V).eulerian_orientation()
    G2 = G.subgraph(V_complement).eulerian_orientation()
    result = DiGraph([G.vertices(), []], multiedges=G.allows_multiple_edges())
    result.add_edges(G.edge_boundary(V))
    result.reverse_edges(result.incoming_edge_iterator(V))
    result.add_edges(G1.edges())
    result.add_edges(G2.edges())
    if show:
        result.show(vertex_colors={'b': V, 'r': V_complement})
    return result
