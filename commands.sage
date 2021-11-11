from time import sleep

def CCSgen(n):
    return (CycleCocycleSystem(G) for G in graphs(n) if G.is_biconnected())

load("CycleCocycleSystem.sage")
G = DiGraph(20)
G.add_cycle(range(12))
G.add_cycle(range(12,18))
G.add_cycle([18,19,20])
G.add_edges([(2*n, n + 12) for n in range(0, 6)])
G.add_edges([(13,20), (17,19), (15,18)])
U = SuperDiGraph(G, G, unori={(i, (i+1) % 12) for i in range(12)})
U.modified_unfurl(range(12))
