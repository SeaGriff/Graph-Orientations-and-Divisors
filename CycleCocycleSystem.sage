""" Builds a cycle-cocycle system on top of a graph and implements associated
methods, in particular those that relate orientations to divisors. """

load("newmethods.sage")

class CycleCocycleSystem(Graph):
    """ Class for a generalized cycle-cocycle system - a graph with
    orientations and related divisorial data. """
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
            self._base_orientation = QuasiDiGraph(self,
                                                  Graph(self).random_orientation())
        else:
            self._base_orientation = QuasiDiGraph(self,
                                                  base_orientation)
        if base_edge is None:
            self._base_edge = self.random_edge()

    """ Basic class functionality. """

    def _repr_(self):
        return ("A graph with a cycle-cocycle reversal system, on {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self):
        Graph(self).show()

    def copy(self):
        return CycleCocycleSystem(self, self._base_orientation, self._base_edge)

    """ Random data """

    def random_orientation(self, unori=0, biori=0):
        U = QuasiDiGraph(self, Graph(self).random_orientation())
        for i in range(unori):
            U.unorient_edge(U.random_oriented_edge())
        for i in range(biori):
            U.biorient_edge(U.random_oriented_edge())
        return U

    """ Set or retrieve attached objects and internal variables """

    def set_base_edge(self, e):
        """ Sets the base edge """
        self._base_edge = e

    def set_base_orientation(self, G):
        """ Sets the base orientation """
        self._base_orientation = G

    def base_orientation(self):
        """ returns the base orientation """
        return self._base_orientation.copy()

    def base_edge(self):
        """ returns the base edge """
        return self._base_edge

    """ Get invariant data """

    def genus(self):
        """ Returns the genus of the graph. """
        return len(self.edges()) - len(self.vertices()) + 1

    def pic(self):
        """ return the picard group of the graph """
        return self._pic

    def orientation_representatives(self):
        """ Returns a list of representatives for the cycle cocycle system """
        return [self.linear_orientation_class(D) for D
                in self._pic.picard_representatives(self.genus() - 1)]

    def theta_char_divisors(self):
        """ Returns a list of all theta characteristic divisors for G. Slow """
        return [D for D in self._pic.picard_representatives(self.genus() -
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
        return QuasiDiGraph(self, G)

    """ Divisorial algorithms """

    def linear_orientation_class(self, div, curU=self.base_orientation()):
    """ This implements an algorithm from section 4 of Backman's 2017
    paper "Riemann-Roch Theory for Graph Orientations."

    Accepts a divisor div and returns an orientation with that div as its
    Chern class, if possible, and an acyclic orientation with Chern class
    dominating div otherwise. """
        assert div.deg() < self.genus(), "Divisor must have degree at most g - 1."
        curD = div - U.chern_class()
        zero = self._pic.zero_div()
        while not curD.is_linearly_equivalent(zero):
            S = curD.div_pos()
            R = (-1*curD).div_pos()
            T = curU.adjacent_to_unori()
            if S != zero and len(curU.unori()) != 0:
                self._linear_orientation_class_case_1(curU, curD, S, T)
            if S != 0 and len(curU.unori()) == 0:
                self._linear_orientation_class_case_2(curU, curD, S, R)
            if S == zero: # Case 3
                incoming = curU.traverser().incoming_edges(R.support())
                if len(incoming) != 0:
                    curU.unorient_edge(incoming[0])
                    curD[incoming[0][1]] += 1
                else:
                    return curU
        return curU

    def _linear_orientation_class_case_1(self, curU, curD, S, T):
        bar_S = curU.reachable_from_vertices(S.support())
        while bar_S.isdisjoint(T):
            curU.reverse_edges(curU.edge_boundary(
                               curU.vertex_complement(S)))
            bar_S = curU.reachable_from_vertices(S.support())
        P = next(curU.traverser().all_paths_iterator(S.support(), T,
                 simple=True, trivial=True))
        curU.reverse_edges(curU._unori.edge_boundary({P[-1]}))
        curU.remove_unorientation(curU._unori.incoming_edges({P[-1]})[0])
        curU.reverse_edges(zip(P, P[1:]))
        curD[P[0]] += -1

    def _linear_orientation_class_case_2(self, curU, curD, S, R):
        bar_S = curU.reachable_from_vertices(S.support())
        if bar_S.isdisjoint(R):
            curU.reverse_edges(curU.edge_boundary(
                               curU.vertex_complement(S)))
            bar_S = curU.reachable_from_vertices(S.support())
        else:
            P = next(curU.traverser().all_paths_iterator(S.support(),
                     R.support(), simple=True))
            curU.reverse_edges(zip(P, P[1:]))
            curD[P[0]] += -1
            curD[P[-1]] += 1

    """ Test associated objects """

    def is_theta_char(self, T):
        """  """
        if isinstance(T, QuasiDiGraph):
            if len(T.biori_set()) != len(T.unori_set()):
                return False
            T = T.chern_class()
        return T.is_linearly_equivalent(div_op(self._pic, T))

class QuasiDiGraph(DiGraph):
    """ Accepts a DiGraph and implements methods to view it as a quasidirected
    graph (a graph where edges may additionally be unoriented, or
    unoriented in both directions).

    Unoriented and bioriented edges are stored in attached graphs,
    _unori and _bi. """


    def __init__(self, ccs, data, bi={}, unori={}):
        DiGraph.__init__(self, data)

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
        return ("A superorientation on a graph with {} vertices and with {} edges".format(
        len(self.vertices()), len(self.edges())))

    def show(self, biori_color='blue', unori_color='red'):
        """ Displays all edges; bioriented edges are by default blue,
        and unoriented edges are by default red. """
        DiGraph(self).show(edge_colors={biori_color: self.biori(),
                                        unori_color: self.unori()})

    def copy(self):
        return QuasiDiGraph(self._ccs, self, self.biori(), self.unori())

    """ Random data """

    def random_oriented_edge(self):
        U = self.copy()
        U._del_unori()
        return U.random_edge()

    def random_unori_edge(self):
        return self._unori.random_edge()

    def random_biori_edge(self):
        return self._biori.random_edge()

    """ Set or retrieve attached objects and variables. """

    def biori(self):
        return self._bi.edges()

    def unori(self):
        return self._unori.edges()


    def traverser(self):
        """ Makes a graph that algorithms can crawl through for paths. """
        U = self.copy()
        U._del_unori()
        U._double_bi()
        return DiGraph(U)

    def ccs(self):
        """ Returns the associated cycle cocycle system. """
        return self._ccs.copy()

    def unorient_edge(self, e, check=True):
        if check:
            assert self.has_edge(e), "Unorienting directed edge not in the orientation."
        self._unori.add_edge(e)

    def biorient_edge(self, e, check=True):
        if check:
            assert self.has_edge(e), "Biorienting directed edge not in the orientation."
        self._bi.add_edge(e)

    def remove_unorientation(self, e):
        self._unori.delete_edge(e)

    def remove_biorientation(self, e):
        self._bi.delete_edge(e)

    def unorient_edges(self, X, check=True):
        for e in X:
            self.unorient_edge(e, check)

    def biorient_edges(self, X, check=True):
        for e in X:
            self.biorient_edge(e, check)

    def set_unori(self, X):
        """ Sets the unoriented edges to be exactly X. """
        self.remove_unorientation(self.unori())
        self.unorient_edges(X)

    def set_biori(self, X):
        """ Sets the bioriented edges to be exactly X. """
        self.remove_biorientation(self.biori())
        self.biorient_edges(X)

    """ Overriding DiGraph methods """

    def reverse_edge(self, e):
        if self._unori.has_edge(e):
            self._unori.reverse_edge(e)
        if self._bi.has_edge(e):
            self._bi.reverse_edge(e)
        DiGraph.reverse_edge(self, e)

    def reverse_edges(self, X):
        for e in X:
            self.reverse_edge(e)

    def reverse(self):
        self._unori.reverse()
        self._bi.reverse()
        DiGraph.reverse(self)

    def sources(self):
        return DiGraph(self.traverser()).sources()

    def reachable_from_vertex(self, q):
        return self.traverser().reachable_from_vertex(q)

    def reachable_from_vertices(self, X):
        return self.traverser().reachable_from_vertices(X)

    """ Test associated objects """

    def is_unoriented(self, e):
        return self._unori.has_edge(e)

    def is_bioriented(self, e):
        return self._bi.has_edge(e)

    def adjacent_to_unori(self):
        return {v for v in self._unori.vertices() if self._unori.degree(v) > 0}

    def adjacent_to_biori(self):
        return {v for v in self._biori.vertices() if self._biori.degree(v) > 0}

    def pivot_into_cut(self, X):
        """ Pivots all possible oriented edges into the boundary of X,
        pointing toward X. """
        edges_it = self.undirected_boundary(X)
        for e in edges_it:
            ind_G = DiGraph(self.traverser()).subgraph(X)
            incoming_at_X_end = ind_G.incoming_edges(e[1])
            if len(incoming_at_X_end) != 0:
                self.edge_pivot(e, incoming_at_X_end[0])

    """ Nice representations """

    def chern_class(self):
        """ returns the Chern class of the orientation. """
        D = self._ccs.pic().all_k_div(-1)
        for e in self.traverser().edges():
            D[e[1]] += 1
        return D

    def q_red(self, q=self.vertices()[0]):
        """ performs the orientation equivalent of passing to the q-reduced
        representative of a divisor """
        self.make_paths(orientation, q)

    def edge_pivot(self, unori_edge, ori_edge):
        """ Performs an edge pivot on an oriented edge (u,v) and an unoriented
        edge adjacent to v. """
        if unori_edge[1] != ori_edge[1]:
            self.reverse_edge(unori_edge)
        self.unorient_edge(ori_edge)
        self.remove_unorientation(unori_edge)

    """ The following implements algorithms from section 4 of Backman's 2017
    paper "Riemann-Roch Theory for Graph Orientations." """

    def dhars(self, early_termination_data=False):
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
        (U, X) = self.dhars(True)
        if X == set(U.vertices()):
            return U
        U.reverse_edges(U.edge_boundary(X))
        return U.unfurl()

    def modified_unfurl(self, S):
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
    """ Checks whether self and another quasiorientation have linearly
    equivalent Chern classes. """
        return self.chern_class().is_linearly_equivalent(U.chern_class())

    def is_theta_char(self):
    """ Checks whether self is a theta character. """
        return self._ccs.is_theta_char(self)

    def undirected_boundary(self, X):
        """ Returns all undirected edges at the boundary of a collection
        of vertices X. """
        self.reverse_edges(self._unori.outgoing_edges(X))
        return self._unori.incoming_edges(X)

    """ Lossy internal methods """

    def _del_unori(self):
        self.delete_edges(self.unori())

    def _double_bi(self):
        self.add_edges(self.biori())
