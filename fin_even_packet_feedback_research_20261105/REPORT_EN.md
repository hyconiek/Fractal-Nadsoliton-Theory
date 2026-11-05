# The even cubic trace budget closes for the finite-packet comparison

EF controls feedback of PC's retarded odd response into the even field.
A duality argument upgrades the odd local-energy estimate to a global
Duhamel energy bound. This allows the broad spatial product with the
packet to be estimated in global L2, without inventing positive moments
or dropping the outgoing harmonic history.

The generated shape response is retained separately and controlled using
its nonresonant scalar carriers and the actual packet envelope. The direct
cubic source is also included. At packet amplitude epsilon and width
R=epsilon^-1, their integrated squared compact-trace bounds are:

| Contribution | Bound |
|---|---|
| Odd-continuum feedback | O(epsilon^2) |
| Shape feedback | O(epsilon^5) |
| Direct cubic term with fixed g4 | O(epsilon^3) |

Thus the combined cubic comparison response meets O(epsilon^2) on the
critical interval. This result concerns TP's specified finite-energy
packet and its retarded response hierarchy. It is not a proof that the
actual nonlinear solution is approximated by that hierarchy or that the
original zero-continuum preparation satisfies the same bound.

Next: a quantified stability/superposition estimate with cubic cross terms,
followed by a justified match to the actual generated continuum. The
separate modal action errors and the actual nonlinear cubic gate remain open.

## Finite numerical comparison

m^2=16/3, S=1; the sum uses g4=1 only as a documented basis example.

| R | dt | dx | Half-domain | Integrated squared total trace |
|---|---|---|---|---|
| 4 | 0.02 | 0.0625 | 64 | 6.29057217469 |
| 4 | 0.01 | 0.0625 | 64 | 6.2899853186 |
| 4 | 0.01 | 0.03125 | 64 | 6.29140063677 |
| 4 | 0.01 | 0.0625 | 96 | 6.2899853186 |
| 8 | 0.02 | 0.0625 | 128 | 69.4163572225 |
| 8 | 0.01 | 0.0625 | 128 | 69.434288355 |

Halving dt reduces each energy/source-work defect by approximately four.
The largest defect in these runs is below 0.00022. Halving dx changes the
compact-observation quadrature by about 0.023%; enlarging the R=4 domain
changes the result only at floating-point precision. The JSON retains the
full three-channel trace Gram matrices, including interference. These six
finite comparisons do not establish a scaling exponent or a nonlinear
approximation theorem.
