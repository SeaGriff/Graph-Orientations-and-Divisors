This project implements algorithms in or connected to
my preprint on the graph Torelli theorem, located at 
https://arxiv.org/abs/2112.07009.

The underlying theory largely deals with generalized cycle
cocycle systems (henceforth just cycle cocycle systems or CCSs).
These were introduced in Spencer Backman's "Riemann-Roch Theory for Graph
Orientations."

The algorithms are of two types:

1) Computations concerning morphisms between CCSs. These are new.

2) Those algorithms appearing in Backman's paper upon which the results I used
depend. I was unable to find other existing implementations of this material.

######

The code here is a foundation: it is not optimized, and various features ought
to be fleshed out. For example, the data structure supports bioriented edges,
but the algorithms have not been adapted for them.

The intent is to produce something others can experiment with or
build upon.

######

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


###


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


bi: An iterable of edges to be initially bioriented.


unori: Likewise for unoriented edges.


###


- OrCycMorphism: Parent is dict.
This is a class for morphisms in the category OrCyc, and in particular
their functorial aspects. Each instance is a dict relating edge labels in a
domain CCS to edge labels in a codomain CCS, encoding a cyclic bijection
(that is, an isomorphism of their matroids).

Produces the signs associated with edges, which are used in functorially
pushing forward orientations and divisors. Computes both forms of pushing
forward.

Checks rigidity, and when rigidity holds, computes the series fixing
automorphism of the main theorem.


Construction:
OrCycMorphism(G, H, f)


G and H: Two CCSes.


f: A dict encoding (by edge labels) a cyclic bijection. If G and H in fact have
such an f, it can produced by the command:
f = Matroid(G).isomorphism(Matroid(H))

####

to test:

At Sage prompt, do

sage: load('example_constructions.sage')
sage: ex1()

It should print

(A graph with a cycle-cocycle reversal system, on 5 vertices and with 7 edges,
 A graph with a cycle-cocycle reversal system, on 5 vertices and with 7 edges,
 A morphism between cycle cocycle systems.)



