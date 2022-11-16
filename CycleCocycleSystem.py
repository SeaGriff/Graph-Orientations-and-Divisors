"""Builds a cycle-cocycle system on top of a graph and implements associated
methods, in particular those that relate orientations to divisors.
"""

from sage.graphs.graph import DiGraph, Graph
from sage.sandpiles.sandpile import Sandpile, SandpileDivisor

from newmethods import *


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
            if len(set(data.edge_labels())) < data.size():
                data = data.autolabel()

        # Build the graph
        Graph.__init__(self, data)

        # Error checking
        if check:
            assert self.is_biconnected(), "Graph is not 2-edge connected."

        # Build associated objects
        self._pic = Sandpile(self)
        if base_orientation is None:
            self._base_orientation = QuasiDiGraph(self.random_orientation(),
                                                  self)
        elif len(set(base_orientation.edge_labels())) < base_orientation.size():
            self._base_orientation = QuasiDiGraph(base_orientation.autolabel(),
                                                  self)
        else:
            self._base_orientation = QuasiDiGraph(base_orientation,
                                                  self)

        if base_edge is None:
            self._base_edge = self.random_edge()
        else:
            self._base_edge = base_edge

    """ Basic class functionality. """

    def _repr_(self):
        return ("A graph with a cycle-cocycle reversal system, on {} vertices and with {} edges".format(
        self.order(), self.size()))

    def show_unoriented(self, edge_labels=True, **kwargs):
        """Show the (unoriented) underlying graph, with edge labels by default."""
        Graph(self).show(edge_labels=edge_labels, **kwargs)

    def show(self, edge_labels=True, **kwargs):
        """Show the base orientation, with edge labels by default."""
        self._base_orientation.show(edge_labels=edge_labels, **kwargs)

    def copy(self):
        """Shallow copy the CCS."""
        return CycleCocycleSystem(self, self._base_orientation, self._base_edge)

    """ Random data """

    def random_orientation(self, unori=0, biori=0):
        """
        Return a random orientation with unori many unoriented edges
        and biori many bioriented edges.
        """
        U = QuasiDiGraph(Graph(self).random_orientation(), self)
        for i in range(unori):
            U.unorient_edge(U.random_oriented_edge())
        for i in range(biori):
            U.biorient_edge(U.random_oriented_edge())
        return U

    """ Set, create, or retrieve attached objects and internal variables """

    def set_base_edge(self, e, label=True):
        """Set the base edge. If label is True, set by label."""
        if label:
            e = self.label_edge_dict()[e]
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
        edge_dict = self.label_edge_dict()
        U.reverse_edges({edge_dict[l] for l in signs.keys() if signs[l] == -1})
        return U

    """ Get invariant data """

    def genus(self):
        """Return the genus of the graph."""
        return self.size() - self.order() + 1

    def pic(self):
        """Return the picard group of the graph."""
        return self._pic

    def orientation_representatives(self):
        """Return a list of representatives for the cycle cocycle system."""
        return [self.orientation_class(D) for D
                in self._pic.picard_representatives(self.genus() - 1)]

    def theta_char_divisors(self):
        """
        Return a list of all theta characteristic divisors for G.
        A theta characteristic is a divisor which is self dual under the
        hodge star operator.
        """
        theta_char = self.sample_theta_char_div()
        return [theta_char + D for D in self._pic.n_torsion(2)]

    def theta_char_orientations(self):
        """Return a list of all theta characteristic orientations for G."""
        return [self.orientation_class(D) for D in
                self.theta_char_divisors()]

    def sample_theta_char_div(self):
        """Return a single theta characteristic divisor. Fast."""
        return self.sample_theta_char_orientation().chern_class()

    def sample_theta_char_orientation(self, show=False):
        """Return a single theta characteristic orientation. Fast."""
        V = Graph(self).eulerian_bipartition()
        return QuasiDiGraph(self._partition_to_ori(V, show), self)

    def _partition_to_ori(V, show=False):
        """
        Accept a nonempty subset V of the vertices. The subset and
        complement are assumed to induce two Eulerian connected subgraphs.
        Return an orientation which is Eulerian on the two subgraphs and
        has a consistently oriented cut between them.
        """
        G = Graph(self)
        if set(V) == set(G.vertices(sort=False)):
            return G.eulerian_orientation()
        V_comp = G.vertex_complement(V)
        G1 = G.subgraph(V).eulerian_orientation()
        G2 = G.subgraph(V_comp).eulerian_orientation()
        result = DiGraph([G.vertices(sort=False), []], multiedges=G.allows_multiple_edges())
        result.add_edges(G.edge_boundary(V))
        result.reverse_edges(result.incoming_edge_iterator(V))
        result.add_edges(G1.edges(sort=False))
        result.add_edges(G2.edges(sort=False))
        if show:
            result.show(vertex_colors={'b': V, 'r': V_comp})
        return result

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

    def orientation_class(self, div, U=None):
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
            return edge_signs(self._base_orientation.edges(sort=False), U.edges(sort=False))
        return edge_signs(self._base_orientation.edges(sort=False), U)


class QuasiDiGraph(DiGraph):
    """
    Implements methods which view a DiGraph as a
    quasidirected graph (a graph where edges may additionally be unoriented, or
    unoriented in both directions).

    Unoriented and bioriented edges are stored in attached graphs,
    _unori and _bi.
    """

    def __init__(self, data, ccs=None, bi=None, unori=None, autolabel=True):

        self._ccs = None

        # Check for distinct edge labels
        if autolabel:
            data = DiGraph(data)
            if len(set(data.edge_labels())) < data.size():
                data = data.autolabel()

        # Construct the quasidigraph
        DiGraph.__init__(self, data)

        if bi is None:
            bi = {}
        if unori is None:
            unori = {}
        self._bi = DiGraph([self.vertices(sort=False), bi],
                           multiedges=self.allows_multiple_edges(),
                           format='vertices_and_edges')
        self._unori = DiGraph([self.vertices(sort=False), unori],
                              multiedges=self.allows_multiple_edges(),
                              format='vertices_and_edges')

        if ccs is None:
            self._ccs = CycleCocycleSystem(Graph(self), base_orientation=self)
        elif isinstance(ccs, CycleCocycleSystem):
            self._ccs = ccs
        else:
            self._ccs = CycleCocycleSystem(ccs)

    """Operator overloading"""

    def __add__(self, D):
        """
        Add a sandpile divisor D (on the right) to self, using the canonical
        pic^0 torsor structure on orientations.
        """
        assert isinstance(D, SandpileDivisor) and D.deg() == 0, "Can only add a degree zero divisor to an orientation."
        assert len(self.unori()) == 0 and len(self.biori()) == 0, "Cannot have unoriented or bioriented edges."
        newD = self.chern_class()
        return self._ccs.orientation_class(D + newD, self)

    def __sub__(self, U):
        """
        Take the difference with U in the canonical pic^0 torsor structure
        on orientations.
        """
        assert len(self.unori()) == 0 and len(self.biori()) == 0, "Cannot have unoriented or bioriented edges."
        assert len(U.unori()) == 0 and len(U.biori()) == 0, "Cannot have unoriented or bioriented edges."
        return self.chern_class() - U.chern_class()

    """Basic class functionality"""

    def __repr__(self):
        """Return a brief description of the class."""
        return ("A quasiorientation on a graph with {} vertices and with {} edges".format(
        self.order(), self.size()))

    def show(self, biori_color='blue', unori_color='red', edge_labels=True,
             **kwargs):
        """
        Display the graph. Bioriented edges are by default blue,
        and unoriented edges are by default red. Edge labels on by default.
        """
        DiGraph(self).show(**kwargs, edge_labels=edge_labels,
                           edge_colors={biori_color: self.biori(),
                           unori_color: self.unori()})

    def copy(self):
        """Shallow copy the quasidigraph."""
        return QuasiDiGraph(self, self._ccs, self.biori(), self.unori())

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
        return self._bi.edges(sort=False)

    def unori(self):
        """Return a list of unoriented edges."""
        return self._unori.edges(sort=False)

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
        """Return the (by default oriented) base edge."""
        return self._ccs.base_edge(oriented)

    def base_label(self):
        """Return the label of the base edge."""
        return self._ccs.base_edge()[2]

    def unorient_edge(self, e, check=True):
        """
        Unorient an edge.
        Buggy behaviour will result from repeatedly unorienting the same edge.
        """
        if check:
            assert self.has_edge(e), "Unorienting directed edge not in the orientation."
            if self._bi.has_edge(e):
                self.remove_biorientation(e)
        self._unori.add_edge(e)

    def biorient_edge(self, e, check=True):
        """
        Biorient an edge.
        Buggy behaviour will result from repeatedly biorienting the same edge.
        """
        if check:
            assert self.has_edge(e), "Biorienting directed edge not in the orientation."
            if self._unori.has_edge(e):
                self.remove_unorientation(e)
        self._bi.add_edge(e)

    def unorient_edges(self, X, check=True):
        """Unorient a collection of edges.
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

    def remove_unorientation(self, e):
        """Give an unoriented edge an (arbitrary) orientation."""
        self._unori.delete_edge(e)

    def remove_biorientation(self, e):
        """Give a bioriented edge an (arbitrary) single orientation."""
        self._bi.delete_edge(e)

    def set_unori(self, X):
        """Set the unoriented edges to be exactly X."""
        for e in self.unori():
            self.remove_unorientation(e)
        self.unorient_edges(X)

    def set_biori(self, X):
        """Set the bioriented edges to be exactly X."""
        for e in self.biori():
            self.remove_biorientation(e)
        self.biorient_edges(X)

    def hodge_star(self):
        """Reverse all edges and swap unorientation and biorientation."""
        self.reverse()
        u, b = self.unori(), self.biori()
        self.set_unori(b)
        self.set_biori(u)

    """ Overriding DiGraph methods """

    def reverse_edge(self, e, multiedges=True):
        """Reverse an edge."""
        if self._unori.has_edge(e):
            self._unori.reverse_edge(e, multiedges)
        if self._bi.has_edge(e):
            self._bi.reverse_edge(e, multiedges)
        if self.has_edge(e):
            DiGraph.reverse_edge(self, e, multiedges=multiedges)
        else:
            DiGraph.reverse_edge(self, (e[1], e[0], e[2]), multiedges=multiedges)

    def reverse_edges(self, X):
        """Reverse a collection of edges."""
        for e in X:
            self.reverse_edge(e)

    def reverse(self):
        """
        Reverses all edges in the graph. Does not affect unoriented and
        bioriented edges.
        """
        for e in self.edges(sort=False):
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
        return {v for v in self._unori.vertices(sort=False) if self._unori.degree(v) > 0}

    def adjacent_to_biori(self):
        """Return a list of all vertices adjacent to a bioriented edge."""
        return {v for v in self._bi.vertices(sort=False) if self._bi.degree(v) > 0}

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
        for e in self.traverser().edges(sort=False):
            D[e[1]] += 1
        return D

    def q_red(self, q=None):
        """
        Perform the orientation equivalent of passing to the q-reduced
        representative of a divisor.
        """
        if q is None:
            q = self.vertices(sort=False)[0]
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
        if X == set(U.vertices(sort=False)):
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
            if X == self.vertices(sort=False):
                return self
            return self._mod_unfurl_it(S, X)
        self.reverse_edges(self.edge_boundary(X))
        if len(self.traverser().incoming_edges(S)) != 0:
            return self
        return self._mod_unfurl_it(S, S)

    """ Test associated objects """

    def equiv(self, U):
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
        return edge_signs(self._ccs.edges(sort=False), self.edges(sort=False))

    def compare_to_self(self, U):
        """
        Return a dict indexed by edge labels in U, comparing the
        edges of U to edges of self. U may be either an orientation
        or a collection of edges. The entries are
        - if edge directions match: 1
        - if edge directions are opposite: -1
        """
        if isinstance(U, DiGraph):
            return edge_signs(self.edges(sort=False), U.edges(sort=False))
        return edge_signs(self.edges(sort=False), U.edges(sort=False))

    """ Lossy internal methods """

    def _del_unori(self):
        self.delete_edges(self.unori())

    def _double_bi(self):
        self.add_edges(self.biori())


class OrCycMorphism(dict):
    """
    A class for morphisms in OrCyc. The inputs are two cycle cocycle systems
    G, H and a morphism f. This f must be formatted as a dict with keys the
    labels of edges in G and entries the labels of edges mapped to.
    The reason for the strange input formatting is that this is
    is the output format SageMath uses when asked for an isomorphism
    Matroid(G) -> Matroid(H) when G and H have labeled edges
    (which CCSes are designed to always have).

    The internal copy of the codomain has its base edge set to match that of
    the domain.

    At construction, determines the correct signs for f.
    Construction is not optimized!

    Has methods for mapping orientations to orientations
    (only full orientations are implemented) and divisors to divisors
    via the Jacobian pushforward.
    """

    def __init__(self, G, H, f):
        """Construct the morphism."""

        dict.__init__(self)
        self.update({l: f[l] for l in f.keys()})

        self._dom = G.copy()
        self._cod = H.copy()
        self._cod.set_base_edge(f[G.base_label()])

        # Construct the map signs, stored in _signs
        cycs = self._orient_cycle_basis(G)
        hcycs = []
        for C in cycs:
            hcycs.append({f[l] for l in C.keys()})
        hcycs = self._orient_cycle_basis(H, hcycs)
        self._signs = {G.base_label(): 1}
        init = next(i for i, C in enumerate(cycs) if G.base_label() in C.keys())
        self._crawl_cycles_comparing(G.base_label(), init,
                                     set({}), cycs, hcycs)


    def _orient_cycle_set(self, CCS, C, align_to=None, is_sorted=False):
        if is_sorted:
            cycle = CCS.compare_to_base_ori(C)
        else:
            edge_set = set({e for e in CCS.edges(sort=False) if e[2] in C})
            T = Graph(CCS).subgraph(edges=edge_set)
            edge_C = T.cycle_basis("edge")[0]
            cycle = CCS.compare_to_base_ori(edge_C)
        if align_to is not None and cycle[align_to] == -1:
            cycle = {l: cycle[l] * -1 for l in cycle.keys()}
        return cycle

    def _crawl_cycles_comparing(self, l, index, checked, cycs,
                                hcycs):
        C = cycs[index]
        fC = hcycs[index]
        pm = self._signs[l] * C[l] * fC[self[l]]
        self._signs.update({l: pm * C[l] * fC[self[l]] for l in C.keys()})
        to_it = []
        for i, newC in enumerate(cycs):
            if len(set(C.keys()).intersection(newC.keys()) - checked) != 0:
                to_it.append(i)
        checked.update(C.keys())
        for i in to_it:
            l = set(C.keys()).intersection(cycs[i].keys()).pop()
            self._crawl_cycles_comparing(l, i, checked, cycs, hcycs)

    def _orient_cycle_basis(self, CCS, cycs=None):
        if cycs is None:
            return [CCS.compare_to_base_ori(C) for C in CCS.cycle_basis(oriented=False)]
        return [self._orient_cycle_set(CCS, C) for C in cycs]

    """Basic class functionality."""

    def __repr__(self):
        return "A morphism between cycle cocycle systems."

    def show(self):
        """Print the underlying dict."""
        print(dict(self))

    def copy(self):
        """Return a shallow copy of self."""
        return OrCycMorphism(self._dom, self._cod, dict(self))

    """Retrieve information about associated objects."""

    def domain(self):
        """Return the domain."""
        return self._dom.copy()

    def codomain(self):
        """Return the codomain."""
        return self._cod.copy()

    def signs(self):
        """Return the dict encoding the morphism's signs."""
        return self._signs.copy()

    def rigidity_divisor(self):
        """Return the rigidity divisor of the morphism."""
        A = self.map(self._dom.base_orientation().chern_class())
        B = self.map(self._dom.base_orientation()).chern_class()
        return A - B

    def is_rigid(self):
        """Return whether the morphism is rigid."""
        return self.rigidity_divisor().is_linearly_equivalent(
                                        self._cod._pic.zero_div())

    """Map associated objects."""

    def map(self, X):
        """Maps a divisor or orientation on the domain."""
        if isinstance(X, QuasiDiGraph):
            return self._map_ori(X)
        if isinstance(X, SandpileDivisor):
            return self._map_div(X)
        raise TypeError("Input must be a quasidigraph or divisor.")

    def _map_div(self, div):
        D = SandpileDivisor(self._dom._pic, div)
        dom_base = self._dom.base_edge()[1]
        cod_base = self._cod.base_edge()[1]
        D[dom_base] -= div.deg()
        Dp = D.div_pos(False)
        Dn = (-1 * D).div_pos(False)
        dom_jac_div = {l: 0 for l in self._dom.edge_labels()}
        for n, p in zip(Dn, Dp):
            P = next(Graph(self._dom).shortest_simple_paths(n, p, report_edges=True,
                                                     labels=True))
            if P[0][0] != n:
                P[0] = (P[0][1], P[0][0], P[0][2])
            i = 1
            while i < len(P):
                if P[i][0] != P[i - 1][1]:
                    P[i] = (P[i][1], P[i][0], P[i][2])
                i += 1
            S = self._dom.compare_to_base_ori(P)
            for l in S.keys():
                dom_jac_div[l] += S[l]
        cod_jac_div = {self[l]: self._signs[l] * dom_jac_div[l]
                       for l in self.keys()}
        edges = self._cod.label_edge_dict()
        result = self._cod._pic.zero_div()
        result[cod_base] += div.deg()
        for e in self._cod.base_orientation().edges(sort=False):
            result[e[1]] += cod_jac_div[e[2]]
            result[e[0]] -= cod_jac_div[e[2]]
        return result

    def _map_ori(self, U):
        """
        Map an orientation U on the domain to an orientation on the codomain.
        """
        if U is None:
            U = self._dom.base_orientation()
        domain_signs = self._dom.compare_to_base_ori(U)
        cod_signs = {self[l]: domain_signs[l] * self._signs[l]
                     for l in self._dom.edge_labels()}
        return self._cod.ori_from_edge_signs(cod_signs)

    def edge_isomorphism(self, check=True):
        """
        Return False if the morphism is not rigid. Otherwise, compute a
        series fixing automorphism of the codomain psi such that postcomposing
        with psi gives an edge isomorphism.
        """
        if check:
            assert self.is_rigid(), "Morphism must be rigid."

        def op(x):  # Swap 0 and 1
            return 1 - x

        G = self._dom.base_orientation()
        H = self._cod.base_orientation()
        H_ref = H.label_edge_dict()
        G_ref = G.label_edge_dict()

        psi = {H.base_label(): H.base_label()}
        v_iso = {G.base_edge()[0]: H.base_edge()[0],
                    G.base_edge()[1]: H.base_edge()[1]}

        checked = {G.base_edge()[2]}
        ends = {G.base_edge()[0], G.base_edge()[1]}

        while checked != set(G.edge_labels()):
            it = {e for e in G.edges(sort=False) if e[0] in ends or e[1] in ends}
            for e in it - {G_ref[s] for s in checked}:
                l = e[2]
                if e[0] in ends:
                    o = 0
                    ends.add(e[1])
                else:
                    o = 1
                    ends.add(e[0])
                s_class = {b[2] for b in H.series_class(H_ref[self[l]])}
                for u in s_class - set(psi.values()):
                    if H_ref[u][o] == v_iso[e[o]]:
                        psi.update({self[l]: u})
                        v_iso.update({e[op(o)]: H_ref[u][op(o)]})
                        break
                    if H_ref[u][op(o)] == v_iso[e[o]]:
                        psi.update({self[l]: u})
                        v_iso.update({e[op(o)]: H_ref[u][o]})
                        break
                checked.add(l)
        return psi
