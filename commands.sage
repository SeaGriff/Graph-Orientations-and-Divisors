def CCSgen(n):
    return (CycleCocycleSystem(G) for G in graphs(n) if G.is_biconnected())

load("orientations.sage")
G = CycleCocycleSystem(graphs.DodecahedralGraph())
sources = []
while len(sources) == 0:
    U = G.random_orientation()
    sources = U.sources()
W = SuperDiGraph(G, U)
# W.unfurl({sources[0]}).show()
