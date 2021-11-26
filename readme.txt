This project implements code concerning processes within or adjacent to
my forthcoming paper on the graph Torelli theorem.

The paper uses the theory of partial orientations put forward in Spencer
Backman's "Riemann-Roch Theory for Graph Orientations." Although Backman's
paper is algorithmically inclined, I was unable to find any implementations
of the material I used. I thought it might be interesting to others to
implement this material (and the bulk of mine) in SageMath.

The implementation is a foundation: it is not optimized and the data structure
supports bioriented edges, but the algorithms have not been adapted for them.

The project consists of two files. The complex code is in CycleCocycleSystem.
A few simple functions and methods which might be independently useful are
broken out into the other file, newmethods.

CycleCocycleSystem consists of three classes:



- CycleCocycleSystem: Parent is Graph.
This a base which organizes and examines objects in OrCyc having the same
underlying graph. The edges in a CCS must always be uniquely labeled,
and the underlying graph must be 2-edge connected. A CCS is not designed
to have edges added or removed after construction.

This implements Backman's algorithm for attempting to lift a divisor
to a partial orientation. Also implements and adapts an algorithm implicit
in a proof of Gallai, in order to efficiently produce divisors which are
up to linear equivalence half the canonical_divisor.

Construction:
CycleCocycleSystem(data, base_orientation=None, base_edge=None,
autolabel=True, check=True)

data: Any data sufficient for SageMath's Graph constructor (in particular,
a Graph or DiGraph).

base_orientation and base_edge: Used throughout. When no input is chosen these
are both assigned randomly.

autolabel: If True (default), checks whether the edges of the underlying graph
are uniquely labeled. If not, labels them with integers 0 to n, in arbitrary
order.

check: If True (default), raises an error if the underlying graph is not
2-edge connected.



- QuasiDiGraph: Parent is DiGraph.
This is a particular instance of a partial orientation sitting over a CCS.
The term "quasidigraph" refers to a digraph in which edges may be unoriented
or bioriented. A QDG instance contains three graphs: itself,
and two graphs unori and biori with the same vertex set as the
main graph. These graphs are used to track which edges have been unoriented
or bioriented.

This implements Backman's oriented Dhar's algorithm, unfurling algorithm,
and modified unfurling algorithm.

Construction:
QuasiDiGraph(data, ccs=None, bi=None, unori=None)

data: Any data sufficient for SageMath's DiGraph constructor (in particular,
a DiGraph).

ccs: The underlying CCS. If none, a CCS is created from the underlying
undirected graph, and given a base orientation with the same edge directions
as your QDG. This can be a more efficient way of producing a CCS.
Note that if the edges a digraph Q are not uniquely labeled, the edges of
the CCS of QuasiDiGraph(Q) will be autolabeled. However, QuasiDiGraph(Q) itself
will not have labeled edges!
For this reason, for such a Q the recommended pattern is:
Q = QuasiDiGraph(Q).ccs().base_orientation()
which produces a correctly labeled QDG.

bi: An iterable of edges to be initially bioriented.

unori: Likewise for unoriented edges.



- OrCycMorphism: Parent is dict.
This is a class for morphisms in the category OrCyc, and in particular
their functorial aspects. Each instance is a dict relating edge labels in a
domain CCS to edge labels in a codomain CCS, encoding a cyclic bijection
(that is, an isomorphism of their matroids).
Produces the signs associated with edges, which are used in functorially
pushing forward orientations and divisors. Does the functorial pushforward.
Calculates rigidity.

Construction:
OrCycMorphism(G, H, f)
G and H: Two CCSes.
f: A dict encoding (by edge labels) a cyclic bijection. If G and H in fact have
such an f, it can produced by the command:
f = Matroid(G).isomorphism(Matroid(H))
