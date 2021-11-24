from time import sleep
load("CycleCocycleSystem.sage")

G = CycleCocycleSystem(graphs.MoserSpindle())
H = CycleCocycleSystem(graphs.MoserSpindle())
G.set_base_edge(0, True)
H.set_base_edge(0, True)
f = Matroid(G).isomorphism(Matroid(G))
G.map_ori(H, f)

"""G = DiGraph(20)
G.add_cycle(range(12))
G.add_cycle(range(12,18))
G.add_cycle([18,19,20])
G.add_edges([(2*n, n + 12) for n in range(0, 6)])
G.add_edges([(13,20), (17,19), (15,18)])
U = QuasiDiGraph(G, G, unori={(i, (i+1) % 12) for i in range(12)})"""
