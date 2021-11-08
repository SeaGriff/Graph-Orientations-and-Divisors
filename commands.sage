def CCSgen(n):
    return (CycleCocycleSystem(G) for G in graphs(n) if G.is_biconnected())

load("orientations.sage")
