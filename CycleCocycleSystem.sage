"""Builds a cycle-cocycle system on top of a graph and implements associated
methods, in particular those that relate orientations to divisors.
"""

#import newmethods
load("newmethods.sage")
from sage.graphs.graph import Graph
from sage.graphs.graph import DiGraph
from sage.sandpiles.sandpile import Sandpile
from sage.sandpiles.sandpile import SandpileDivisor


class CycleCocycleSystem(Graph):
    """
    Class for a generalized cycle-cocycle system - a graph with
    orientations and related divisorial data.
    """

    def __init__(self, data, base_orientation=None, base_edge=None,
                 autolabel=True, check=True):
        """
        Construct the graph and associated data.
        The data parameter can be any input sufficient for Sage to build
        a graph with.
        """

        # Check for distinct edge labels
        if autolabel:
            data = Graph(data)
            if len(set(data.edge_labels())) < len(data.edges()):
                data = data.autolabel()

        # Build the graph
        Graph.__init__(self, data)

        # Error checking
        if check:
            assert self.is_biconnected(), "Graph is not 2-edge connected."

        # Build associated objects
        self._pic = Sandpile(self)
        if base_orientation is None:
            self._base_orientation = QuasiDiGraph(self,
                                              self.random_orientation())
        else:
            self._base_orientation = QuasiDiGraph(self,
                                                  base_orientation)
        if base_edge is None:
            self._base_edge = self.random_edge()
        else:
            self._base_edge = base_edge

    """ Basic class functionality. """

    def _repr_(self):
        return ("A graph with a cycle-cocycle reversal system, on {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self, **kwargs):
        """Show the (unoriented) underlying graph."""
        Graph(self, **kwargs).show()

    def copy(self):
        """Shallow copy the CCS."""
        return CycleCocycleSystem(self, self._base_orientation, self._base_edge)

    """ Random data """

    def random_orientation(self, unori=0, biori=0):
        """
        Return a random orientation with unori many unoriented edges
        and biori many bioriented edges.
        """
        U = QuasiDiGraph(self, Graph(self).random_orientation())
        for i in range(unori):
            U.unorient_edge(U.random_oriented_edge())
        for i in range(biori):
            U.biorient_edge(U.random_oriented_edge())
        return U

    """ Set, create, or retrieve attached objects and internal variables """

    def set_base_edge(self, e):
        """Set the base edge."""
        self._base_edge = e

    def set_base_orientation(self, G):
        """Set the base orientation."""
        self._base_orientation = QuasiDiGraph(G)

    def base_orientation(self):
        """Return the base orientation."""
        return self._base_orientation.copy()

    def base_edge(self, oriented=True):
        """Return the base edge."""
        if oriented:
            if self._base_orientation.has_edge(self._base_edge):
                return self._base_edge
            else:
                return (self._base_edge[1], self._base_edge[0],
                        self._base_edge[2])
        return self._base_edge

    def base_label(self):
        """Return the label of the base edge."""
        return self._base_edge[2]

    def ori_from_edge_signs(self, signs):
        """
        Accept a dict with keys which are edge labels of self and entries +/-1.
        Return an orientation which differs from the base ori according to the
        signs of entries.
        """
        U = self.base_orientation()
        edge_dict = {e[2]: e for e in self.edges()}
        U.reverse_edges({edge_dict[l] for l in signs.keys() if signs[l] == -1})
        return U

    """ Get invariant data """

    def genus(self):
        """Return the genus of the graph."""
        return len(self.edges()) - len(self.vertices()) + 1

    def pic(self):
        """Return the picard group of the graph."""
        return self._pic

    def orientation_representatives(self):
        """Return a list of representatives for the cycle cocycle system."""
        return [self.linear_orientation_class(D) for D
                in self._pic.picard_representatives(self.genus() - 1)]

    def theta_char_divisors(self):
        """Return a list of all theta characteristic divisors for G."""
        theta_char = self.sample_theta_char_div()
        return [theta_char + D for D in self._pic.n_torsion(2)]

    def theta_char_orientations(self):
        """Return a list of all theta characteristic orientations for G."""
        return [self.linear_orientation_class(D) for D in
                self.theta_char_divisors()]

    def sample_theta_char_div(self):
        """Return a single theta characteristic divisor. Fast."""
        return self.sample_theta_char_orientation().chern_class()

    def sample_theta_char_orientation(self, show=False):
        """Return a single theta characteristic orientation. Fast."""
        V = Graph(self).eulerian_bipartition()
        G = _partition_to_theta_char_orientation(Graph(self), V, show)
        return QuasiDiGraph(self, G)

    def cycle_basis(self, oriented=True):
        """
        Modify the cycle basis method for cycle cocycle systems.
        First, the method now always returns edge formatted cycles.
        Second, add a parameter "oriented".
        When this parameter is True, return the cycle basis
        for the base orientation.
        """
        return self._base_orientation.cycle_basis(oriented)

    """ Divisorial algorithms """

    def linear_orientation_class(self, div, U=None):
        """
        Implements an algorithm from section 4 of Backman's 2017
        paper "Riemann-Roch Theory for Graph Orientations."

        Accepts a divisor div and returns an orientation with that div as its
        Chern class, if possible, and an acyclic orientation with Chern class
        dominating div otherwise.
        """
        assert div.deg() < self.genus(), "Divisor must have degree at most g - 1."
        if U is None:
            U = self.base_orientation()
        else:
            U = U.copy()
        D = div - U.chern_class()
        zero = self._pic.zero_div()
        while not D.is_linearly_equivalent(zero):
            S, R, T = D.div_pos(), (-1*D).div_pos(), U.adjacent_to_unori()
            if S != zero and len(T) != 0: # Case 1
                bar_S = U.reachable_from_vertices(S.support())
                while bar_S.isdisjoint(T):
                    U.reverse_edges(U.edge_boundary(
                                    U.vertex_complement(S.support())))
                    bar_S = U.reachable_from_vertices(S.support())
                P = next(U.traverser().all_paths_iterator(S.support(), T,
                         simple=True, trivial=True))
                U.reverse_edges(U._unori.edge_boundary({P[-1]}))
                U.remove_unorientation(U._unori.incoming_edges({P[-1]})[0])
                U.reverse_edges(zip(P, P[1:]))
                D[P[0]] += -1
            if S != 0 and len(T) == 0: # Case 2
                bar_S = U.reachable_from_vertices(S.support())
                if bar_S.isdisjoint(R.support()):
                    U.reverse_edges(U.edge_boundary(
                                       U.vertex_complement(bar_S)))
                    bar_S = U.reachable_from_vertices(S.support())
                else:
                    P = next(U.traverser().all_paths_iterator(S.support(),
                             R.support(), simple=True))
                    U.reverse_edges(zip(P, P[1:]))
                    D[P[0]] += -1
                    D[P[-1]] += 1
            if S == zero:  # Case 3
                incoming = U.traverser().incoming_edges(R.support())
                if len(incoming) != 0:
                    U.unorient_edge(incoming[0])
                    D[incoming[0][1]] += 1
                else:
                    return U
        return U


    """ Test associated objects """

    def is_theta_char(self, T):
        """Check whether an orientation is a theta characteristic."""
        if isinstance(T, QuasiDiGraph):
            if len(T.biori_set()) != len(T.unori_set()):
                return False
            T = T.chern_class()
        return T.is_linearly_equivalent(T.hodge_star())

    """Functionality for maps between CCSs"""

    def compare_to_base_ori(self, U):
        """
        Accept U, which may be either an orientation or a collection of edges
        of self.
        Return a dict by labels of edges of the base orientation, comparing the
        edges of U to edges of the base orientation. The entries are
        - if edge directions match: 1
        - if edge directions opposite: -1
        """
        if isinstance(U, DiGraph):
            return edge_signs(self._base_orientation.edges(), U.edges())
        return edge_signs(self._base_orientation.edges(), U)

    def orient_cyclic_bijection(self, H, f):
        """
        Accept a cycle cocycle system and a cyclic bijection f. This f must be
        formatted as a dict with keys the labels of edges in self and entries
        the labels of edges mapped to.
        Return f except the entries are tuples (l, s), where s = (+/-)1
        indicates whether f reverses the edge in question.
        The reason for the strange input formatting is that this is
        is the output format SageMath uses when asked for an isomorphism
        Matroid(self) -> Matroid(H) when self and H have labeled edges.

        This function is not optimized.
        """
        H_ref = H.label_edge_dict()
        cycs = self._orient_cycle_basis()
        hcycs = []
        for C in cycs:
            hcycs.append(H.compare_to_base_ori([H_ref[f[l]] for l in C.keys()]))
        hcycs = H._orient_cycle_basis(hcycs)
        result = {self.base_label(): (f[self.base_label()], 1)}
        init = next(i for i, C in enumerate(cycs) if self.base_label() in C.keys())
        self._crawl_cycles_comparing(self.base_label(), init,
                                     set({self.base_label()}), cycs, hcycs,
                                     f, result)
        return result

        def _orient_cycle_basis(self, init_cycs=None):
            if init_cycs is None:
                return [self.compare_to_base_ori(C) for C in self.cycle_basis(oriented=False)]
            else:
                return [self.compare_to_base_ori(C) for C  in init_cycs]

    """def _orient_cycle_basis(self, init_cycs=None):
        if init_cycs is None:
            cycs = [self.compare_to_base_ori(C) for C in self.cycle_basis(oriented=False)]
        else:
            cycs = init_cycs
        initialC = next(C for C in cycs if self.base_label() in C.keys())
        self._crawl_cycles_aligning(self.base_label(), 1, initialC,
                              set({self.base_label()}), cycs)
        return cycs

    def _crawl_cycles_aligning(self, known_l, known_val, C, checked, cycs):
        if C[known_l] != known_val:
            i = cycs.index(C)
            for l in C.keys():
                C[l] *= -1
            cycs[i] = C
        to_it = []
        for newC in cycs:
            if len(set(C.keys()).intersection(newC.keys()) - checked) != 0:
                to_it.append(newC)
        checked.update(C.keys())
        for newC in to_it:
            l = set(C.keys()).intersection(newC.keys()).pop()
            self._crawl_cycles_aligning(l, C[l], newC, checked, cycs)"""

    def _crawl_cycles_comparing(self, l, index, checked, cycs,
                                hcycs, f, result):
        C = cycs[index]
        fC = hcycs[index]
        pm = result[l][1] * C[l] * fC[f[l]]
        result.update({l: (f[l], pm * C[l] * fC[l]) for l in C.keys()})
        to_it = []
        for i, newC in enumerate(cycs):
            if len(set(C.keys()).intersection(newC.keys()) - checked) != 0:
                to_it.append(i)
        checked.update(C.keys())
        for i in to_it:
            l = set(C.keys()).intersection(cycs[i].keys()).pop()
            self._crawl_cycles_comparing(l, i, checked, cycs, hcycs, f, result)

class QuasiDiGraph(DiGraph):
    """
    Implements methods which view a DiGraph as a
    quasidirected graph (a graph where edges may additionally be unoriented, or
    unoriented in both directions).

    Unoriented and bioriented edges are stored in attached graphs,
    _unori and _bi.
    """

    def __init__(self, ccs, data, bi=None, unori=None):
        """Construct the quasidigraph."""
        DiGraph.__init__(self, data)

        if bi is None:
            bi = {}
        if unori is None:
            unori = {}
        self._bi = DiGraph([self.vertices(), bi],
                           multiedges=self.allows_multiple_edges(),
                           format='vertices_and_edges')
        self._unori = DiGraph([self.vertices(), unori],
                              multiedges=self.allows_multiple_edges(),
                              format='vertices_and_edges')

        if isinstance(ccs, CycleCocycleSystem):
            self._ccs = ccs
        else:
            self._ccs = CycleCocycleSystem(ccs)

    """ Basic class functionality """

    def __repr__(self):
        """Return a brief description of the class."""
        return ("A quasiorientation on a graph with {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self, biori_color='blue', unori_color='red', **kwargs):
        """
        Display the graph. Bioriented edges are by default blue,
        and unoriented edges are by default red.
        """
        DiGraph(self).show(**kwargs, edge_colors={biori_color: self.biori(),
                                                  unori_color: self.unori()})

    def copy(self):
        """Shallow copy the quasidigraph."""
        return QuasiDiGraph(self._ccs, self, self.biori(), self.unori())

    """Graph invariants"""

    def cycle_basis(self, oriented=True):
        """
        Modify the cycle basis method for quasidigraphs.
        First, the method now always returns edge formatted cycles.
        Second, add a parameter "oriented".
        When this parameter is True, return the cycle basis
        for self.
        """
        if oriented:
            return DiGraph(self).cycle_basis("edge")
        return Graph(self).cycle_basis("edge")

    """ Random data """

    def random_oriented_edge(self):
        """Return a random oriented (but not bioriented) edge."""
        U = self.copy()
        U._del_unori()
        return U.random_edge()

    def random_unori_edge(self):
        """Return a random unoriented edge."""
        return self._unori.random_edge()

    def random_biori_edge(self):
        """Return a random bioriented edge."""
        return self._bi.random_edge()

    """Set or retrieve attached objects and variables."""

    def biori(self):
        """Return a list of bioriented edges."""
        return self._bi.edges()

    def unori(self):
        """Return a list of unoriented edges."""
        return self._unori.edges()

    def traverser(self):
        """
        Construct a graph that algorithms can crawl through for paths,
        by deleting unoriented edges and doubling bioriented edges.
        """
        U = self.copy()
        U._del_unori()
        U._double_bi()
        return DiGraph(U)

    def ccs(self):
        """Return the associated cycle cocycle system."""
        return self._ccs.copy()

    def base_orientation(self):
        """
        Return the base orientation associated with the underlying cycle
        cocycle system.
        """
        return self._ccs.base_orientation()

    def base_edge(self, oriented=False):
        """Return the (oriented) base edge."""
        return self._ccs.base_edge(oriented)

    def unorient_edge(self, e, check=True):
        """
        Unorient an edge. Buggy behaviour will result from repeatedly
        unorienting the same edge.
        """
        if check:
            assert self.has_edge(e), "Unorienting directed edge not in the orientation."
            if self._bi.has_edge(e):
                self.remove_biorientation(e)
        self._unori.add_edge(e)

    def biorient_edge(self, e, check=True):
        """
        Biorient an edge. Buggy behaviour will result from repeatedly
        biorienting the same edge.
        """
        if check:
            assert self.has_edge(e), "Biorienting directed edge not in the orientation."
            if self._unori.has_edge(e):
                self.remove_unorientation(e)
        self._bi.add_edge(e)

    def remove_unorientation(self, e):
        """Give an unoriented edge an (arbitrary) orientation."""
        self._unori.delete_edge(e)

    def remove_biorientation(self, e):
        """Give a bioriented edge an (arbitrary) single orientation."""
        self._bi.delete_edge(e)

    def unorient_edges(self, X, check=True):
        """U
        norient a collection of edges.
        Buggy behaviour will result from repeatedly unorienting the same edge.
        """
        for e in X:
            self.unorient_edge(e, check)

    def biorient_edges(self, X, check=True):
        """
        Biorient a collection of edges.
        Buggy behaviour will result from repeatedly biorienting the same edge.
        """
        for e in X:
            self.biorient_edge(e, check)

    def set_unori(self, X):
        """Set the unoriented edges to be exactly X."""
        self.remove_unorientation(self.unori())
        self.unorient_edges(X)

    def set_biori(self, X):
        """Set the bioriented edges to be exactly X."""
        self.remove_biorientation(self.biori())
        self.biorient_edges(X)

    def hodge_star(self):
        """Reverse all edges and swap unorientation and biorientation."""
        self.reverse()
        u, b = self.unori(), self.biori()
        self.set_unori(b)
        self.set_biori(u)

    """ Overriding DiGraph methods """

    def reverse_edge(self, u, v=None, label=None, multiedges=True):
        """Reverse an edge."""
        if self._unori.has_edge(u, v, label):
            self._unori.reverse_edge(u, v, label, multiedges)
        if self._bi.has_edge(u, v, label):
            self._bi.reverse_edge(u, v, label, multiedges)
        DiGraph.reverse_edge(self, u, v, label, multiedges)

    def reverse_edges(self, X):
        """Reverse a collection of edges."""
        for e in X:
            self.reverse_edge(e)

    def reverse(self):
        """
        Reverses all edges in the graph. Does not affect unoriented and
        bioriented edges.
        """
        for e in self.edges():
            self.reverse_edge(e)

    def sources(self):
        """Return a list of sources in the graph."""
        return DiGraph(self.traverser()).sources()

    def reachable_from_vertex(self, q):
        """Return a list of vertices reachable by oriented paths from q."""
        return self.traverser().reachable_from_vertex(q)

    def reachable_from_vertices(self, X):
        """
        Return a list of vertices reachable by oriented paths from
        vertices in X.
        """
        return self.traverser().reachable_from_vertices(X)

    """ Test associated objects """

    def is_unoriented(self, e):
        """Return whether an edge is unoriented."""
        return self._unori.has_edge(e)

    def is_bioriented(self, e):
        """Return whether an edge is bioriented."""
        return self._bi.has_edge(e)

    def adjacent_to_unori(self):
        """Return a list of all vertices adjacent to an unoriented edge."""
        return {v for v in self._unori.vertices() if self._unori.degree(v) > 0}

    def adjacent_to_biori(self):
        """Return a list of all vertices adjacent to a bioriented edge."""
        return {v for v in self._bi.vertices() if self._bi.degree(v) > 0}

    def pivot_into_cut(self, X):
        """
        Pivots all possible oriented edges into the boundary of X,
        pointing toward X.
        """
        edges_it = self.undirected_boundary(X)
        for e in edges_it:
            ind_G = DiGraph(self.traverser()).subgraph(X)
            incoming_at_X_end = ind_G.incoming_edges(e[1])
            if len(incoming_at_X_end) != 0:
                self.edge_pivot(e, incoming_at_X_end[0])

    """ Nice representations """

    def chern_class(self):
        """Return the Chern class of the orientation."""
        D = self._ccs.pic().all_k_div(-1)
        for e in self.traverser().edges():
            D[e[1]] += 1
        return D

    def q_red(self, q=None):
        """
        Perform the orientation equivalent of passing to the q-reduced
        representative of a divisor.
        """
        if q is None:
            q = self.vertices()[0]
        self.make_paths(q)

    def edge_pivot(self, unori_edge, ori_edge):
        """
        Perform an edge pivot on an oriented edge (u,v) and an unoriented
        edge adjacent to v.
        """
        if unori_edge[1] != ori_edge[1]:
            self.reverse_edge(unori_edge)
        self.unorient_edge(ori_edge)
        self.remove_unorientation(unori_edge)

    """
    The following implements algorithms from section 4 of Backman's 2017
    paper "Riemann-Roch Theory for Graph Orientations."

    None of these check for correctness of input.
    """

    def dhars(self, early_termination_data=False):
        """
        Accept self with a directed cycle and a source, and
        returns an equivalent orientation which is either acyclic or
        certifies that every equivalent orientation
        contains a directed cycle.
        """
        U = self.copy()
        return U._dhars_it(set(self.sources()), early_termination_data)

    def _dhars_it(self, X, early_termination_data):
        X_comp = self.vertex_complement(X)
        self.pivot_into_cut(X_comp)
        ind_G = self.traverser().subgraph(X_comp)
        v_boundary = self.vertex_boundary(X)
        to_add = {v for v in v_boundary if len(ind_G.incoming_edges(v)) == 0}
        if len(to_add) == 0:
            if early_termination_data:
                return (self, X)
            return self
        return self._dhars_it(X.union(to_add), early_termination_data)

    def unfurl(self):
        """
        Accept self with a directed cycle and a source, and
        return an equivalent orientation which is either acyclic or
        sourceless.
        """
        (U, X) = self.dhars(True)
        if X == set(U.vertices()):
            return U
        U.reverse_edges(U.edge_boundary(X))
        return U.unfurl()

    def modified_unfurl(self, S):
        """
        Accept a collection of sources S inducing a connected
        subgraph, and return
        i) an equivalent orientation with an edge oriented
        into S
        ii) an acyclic orientation in which all the vertices in S are still
        sources. This certifies that in every equivalent orientation,
        the vertices of S are still sources.
        """
        U = self.copy()
        return U._mod_unfurl_it(set(S), set(S))

    def _mod_unfurl_it(self, S, X):
        X_comp = self.vertex_complement(X)
        self.pivot_into_cut(X_comp)
        unori_in_cut = self._unori.edge_boundary(X)
        if len(unori_in_cut) != 0:
            X.add(unori_in_cut[0][1])
            if X == self.vertices():
                return self
            return self._mod_unfurl_it(S, X)
        self.reverse_edges(self.edge_boundary(X))
        if len(self.traverser().incoming_edges(S)) != 0:
            return self
        return self._mod_unfurl_it(S, S)

    """ Test associated objects """

    def is_equivalent(self, U):
        """
        Check whether self and another quasiorientation have linearly
        equivalent Chern classes.
        """
        return self.chern_class().is_linearly_equivalent(U.chern_class())

    def is_theta_char(self):
        """Check whether self is a theta character."""
        return self._ccs.is_theta_char(self)

    def undirected_boundary(self, X):
        """
        Return all undirected edges at the boundary of a collection
        of vertices X.
        """
        self.reverse_edges(self._unori.outgoing_edges(X))
        return self._unori.incoming_edges(X)

    def compare_to_base_ori(self):
        """
        Return a dict indexed by edge labels, comparing the
        edges of self to edges of the base orientation. The entries are
        - if edge directions match: 1
        - if edge directions opposite: -1
        """
        return edge_signs(self._ccs.edges(), self.edges())

    def compare_to_self(self, U):
        """
        Return a dict indexed by edge labels in U, comparing the
        edges of U to edges of self. U may be either an orientation
        of a collection of edges. The entries are
        - if edge directions match: 1
        - if edge directions opposite: -1
        """
        if isinstance(U, DiGraph):
            return edge_signs(self.edges(), U.edges())
        return edge_signs(self.edges(), U.edges())

    """ Lossy internal methods """

    def _del_unori(self):
        self.delete_edges(self.unori())

    def _double_bi(self):
        self.add_edges(self.biori())


def _partition_to_theta_char_orientation(G, V, show=False):
    """
    Accept a graph with a nonempty subset of the vertices. The subset and
    complement are assumed to induce two Eulerian subgraphs.
    Return an orientation which is Eulerian on the two subgraphs and
    has a consistently oriented cut between them.
    """
    if set(V) == set(G.vertices()):
        return G.eulerian_orientation()
    V_comp = G.vertex_complement(V)
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
