"""Here are some example constructions drawn from my paper and Backman's."""


def ex1():
    """
    An example of a rigid morphism. Usage example:

    G, H, phi = ex1()
    G.show()
    H.show()
    phi.show()
    phi.signs()
    phi.is_rigid()
    phi.edge_isomorphism()
    """

    G = QuasiDiGraph({1: {2: [3], 5: [2]}, 2: {3: [4,5]}, 3: {4: [6]}, 5: {1: [1], 4: [7]}}).ccs()
    H = QuasiDiGraph({1: {2: [7], 5: [5]}, 2: {3: [2]}, 3: {2: [1], 4: [6]}, 4: {5: [3]}, 5: {1: [4]}}).ccs()
    G.set_base_edge(1)
    H.set_base_edge(1)
    phi = OrCycMorphism(G, H, {l: l for l in G.edge_labels()})
    return (G, H, phi)


def ex2():
    """
    An example of a nonrigid morphism. Usage example:

    J, K, rho = ex2()
    J.show()
    K.show()
    rho.show()
    rho.signs()
    rho.is_rigid()
    """

    J = QuasiDiGraph({1: {2:[1], 4:[5,6]}, 2: {3:[2,3]}, 3: {4:[4]}}).ccs()
    K = QuasiDiGraph({1: {2:[1], 4:[6]}, 2: {3:[4]}, 3: {4:[2]}, 4: {3:[3], 1:[5]}}).ccs()
    J.set_base_edge(1)
    K.set_base_edge(1)
    rho = OrCycMorphism(J, K, {l: l for l in J.edge_labels()})
    return (J, K, rho)


def ex3():
    """
    The partially oriented graph that Backman unfurls in his paper.
    Usage example:

    U = ex3()
    U.show(edge_labels=False)
    U.unfurl()
    U.show(edge_labels=False)
    """

    U = DiGraph(20)
    U.add_cycle(range(12))
    U.add_cycle(range(12,18))
    U.add_cycle([18,19,20])
    U.add_edges([(2*n, n + 12) for n in range(0, 6)])
    U.add_edges([(13,20), (17,19), (15,18)])
    U = QuasiDiGraph(U)
    U.set_unori(DiGraph(U).subgraph(range(12)).edges())
    return U
