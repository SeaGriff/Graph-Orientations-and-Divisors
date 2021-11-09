""" Builds a cycle-cocycle system on top of a graph and implements associated
methods, in particular those that relate orientations to divisors. """

import itertools
from sage.graphs.graph import Graph
from copy import copy

class SuperDiGraph(DiGraph):
    """ Accepts a DiGraph and implements methods to view it as a super
    directed graph (a graph where edges may additionally be unoriented, or
    unoriented in both directions) """
    def __init__(self, CCS, data, bi={}, unori={}):
        DiGraph.__init__(self, data)

        self._bi = DiGraph(self.vertices(),
                           multiedges=self.allows_multiple_edges())
        self._unori = DiGraph(self.vertices(),
                           multiedges=self.allows_multiple_edges())
        self._bi.add_edges(bi)
        self._unori.add_edges(unori)

        if isinstance(CCS, CycleCocycleSystem):
            self._CCS = CCS
        else:
            self._CCS = CycleCocycleSystem(CCS)

    def __repr__(self):
        return ("A superorientation on a graph with {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self):
        DiGraph(self).show(edge_colors={'blue': self.biori(),
                                        'red': self.unori()})

    def copy(self):
        return SuperDiGraph(self._CCS, self, self.biori(), self.unori())

    def biori(self):
        return self._bi.edges()

    def unori(self):
        return self._unori.edges()

    def CCS(self):
        return self._CCS.copy()

    def unorient_edge(self, e):
        self._unori.add_edge(e)

    def biorient_edge(self, e):
        self._bi.add_edge(e)

    def remove_unorientation(self, e):
        self._unori.remove(e)

    def remove_biorientation(self, e):
        self._bi.remove(e)

    def unorient_edges(self, X):
        for e in X:
            self.unorient_edge(e)

    def biorient_edges(self, X):
        for e in X:
            self.biorient_edge(e)

    def set_unori(self, X):
        self._unori = list(X)

    def set_biori(self, X):
        self._bi = list(X)

    def reverse_edge(self, e):
        e = (e[0], e[1], e[2])
        if e in self._unori:
            self.remove_unorientation(e)
            self.unorient_edge((e[1], e[0], e[2]))
        if e in self._bi:
            self.remove_biorientation(e)
            self.biorient_edge((e[1], e[0], e[2]))
        DiGraph.reverse_edge(self, e) # :< :< :<

    def reverse_edges(self, X):
        for e in X:
            self.reverse_edge(e)

    def op(self):
        return SuperDiGraph(self._CCS, DiGraph(self).reverse(),
                            [(e[1],e[0]) for e in self._unori],
                            [(e[1],e[0]) for e in self._bi])

    def is_equivalent(self, U):
        return self.chern_class().is_linearly_equivalent(U.chern_class())

    def sources(self):
        return DiGraph(self._make_algo_graph()).sources()

    def _pivot_toward(self, X):
        edges_it = self._undirected_boundary(X)
        pre_G = self.copy()
        pre_G._del_unori()
        ind_G = DiGraph(pre_G).subgraph(X)
        for e in edges_it:
            if e[1] not in X:
                self.reverse_edge(e)
            incoming_at_X_end = {l for l in ind_G.incoming_edges(e[1])}
            if len(incoming_at_X_end) != 0:
                self._edge_pivot(e, next(iter(incoming_at_X_end)))

    def dhars(self):
        U = self.copy()
        sources = set(self.sources())
        if self.is_directed_acyclic() or len(sources) == 0:
            return U
        return U.dhars_it(sources)

    def dhars_it(self, X):
        X_complement = set(self.vertices()) - X
        self._pivot_toward(X_complement)
        ind_G = self._make_algo_graph().subgraph(X_complement)
        v_boundary = self.vertex_boundary(X)
        to_add = {v for v in v_boundary if len(ind_G.incoming_edges(v)) == 0}
        if len(to_add) == 0:
            self.show()
            return self
        return self.dhars_it(X.union(to_add))

    # none of the unfurling stuff is done

    def unfurl(self, source):
        U = self.copy()
        sources = set(self.sources())
        if self.is_directed_acyclic() or len(sources) == 0:
            return U

    def modified_unfurl(self, S):
        U = self.copy()
        return self._mod_unfurl_it(set(S), set(S), U)

    def _mod_unfurl_it(self, S, X):
        todo = self._undirected_boundary(X)
        if len(todo) != 0:
            X.add(todo[0][1])
            return self._mod_unfurl_it(S, X, U)
        self.reverse_edges(self.outgoing_edges(X))
        if len(set(self.incoming_edges(S)) - set(self.unori())) != 0:
            return self
        return self._mod_unfurl_it(S, S, self)

    def _edge_pivot(self, unori_edge, ori_edge):
        self.unorient_edge(ori_edge)
        self.remove_unorientation(unori_edge)
        if ori_edge[1] != unori_edge[1]:
            self.reverse_edge(unori_edge)

    def _undirected_boundary(self, X):
        return [e for e in self._unori if e[0] in X and e[1] not in X]

    def _del_unori(self):
        self.delete_edges(self._unori)

    def _double_bi(self):
        self.add_edges({(e[1], e[0]) for e in self._bi})

    def _make_algo_graph(self):
        U = self.copy()
        U._del_unori()
        U._double_bi()
        return DiGraph(U)

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
            self._base_edge = self.random_edge()

    def _repr_(self):
        return ("A graph with a cycle-cocycle reversal system, on {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self):
        Graph(self).show()

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
                        eulerian_bipartition(Graph(self)), show)
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

def double_directed_edges(X):
    """ Accepts an iterable of ordered pairs and returns a set with all
    the pairs and all their reversed versions """
    return set(X).union({(e[1], e[0]) for e in X})

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
