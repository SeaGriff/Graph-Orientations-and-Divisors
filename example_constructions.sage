"""Here are some example constructions drawn from my paper and Backman's:"""

# The first example in my paper:
J = QuasiDiGraph({1: {2:[1], 4:[5,6]}, 2: {3:[2,3]}, 3: {4:[4]}}).ccs()
K = QuasiDiGraph({1: {2:[1], 4:[6]}, 2: {3:[4]}, 3: {4:[2]}, 4: {3:[3], 1:[5]}}).ccs()
J.set_base_edge(1, True)
rho = OrCycMorphism(J, K, {l: l for l in J.edge_labels()})

# The second example in my paper:
G = QuasiDiGraph({1: {2: [3], 5: [2]}, 2: {3: [4,5]}, 3: {4: [6]}, 5: {1: [1], 4: [7]}}).ccs()
H = QuasiDiGraph({1: {2: [7], 5: [5]}, 2: {3: [2]}, 3: {2: [1], 4: [6]}, 4: {5: [3]}, 5: {1: [4]}}).ccs()
G.set_base_edge(1, True)
phi = OrCycMorphism(G, H, {l: l for l in G.edge_labels()})

# The graph that Backman unfurls in his paper
U = DiGraph(20)
U.add_cycle(range(12))
U.add_cycle(range(12,18))
U.add_cycle([18,19,20])
U.add_edges([(2*n, n + 12) for n in range(0, 6)])
U.add_edges([(13,20), (17,19), (15,18)])
G = CycleCocycleSystem(U, base_orientation=U)
U = G.base_orientation()
U.set_unori(DiGraph(U).subgraph(range(12)).edges())
