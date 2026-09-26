# TWO-SCALE-METASTABLE-RESISTANCE-GEOMETRY-116
## Fast four-cycles and a slow three-state quotient form an exact asymptotic two-scale network

Date: 2026-09-26

Status:
- exact graph reduction in the hierarchy k3 >> k4 >> remaining rates;
- rate exponents inherited from reports 110-112;
- no continuum or physical dimension is inferred.

## 1. Fast graph

The cheapest d=3 transitions split the twelve minima into

    C0={0,3,6,9},
    C1={1,4,7,10},
    C2={2,5,8,11}.

Each C_a is a four-cycle.

At the fast scale each cluster internally equilibrates while the other two
clusters are effectively disconnected.

## 2. Slow quotient

Each microscopic state has:
- one d=4 edge into one other residue class;
- one d=4 edge into the remaining residue class.

After the fast C4 equilibration, exact lumping gives the three-state generator

        [-2k4,  k4,  k4]
    Q = [ k4, -2k4,  k4]
        [ k4,  k4, -2k4].

Its nonzero decay eigenvalues are

    boxed:
    3 k4, 3 k4.

This is the same slow pair already identified as Fourier modes m=4,8.

## 3. Coarse conductance

The quotient stationary distribution is uniform:

    pi_cluster=1/3.

Hence each quotient edge has reversible conductance

    C_cluster
      = k4/3.

The three coarse sectors therefore form an equilateral conductance triangle.

Their effective electrical resistance is

    R_eff(cluster a,cluster b)
      = 2/k4.

So the slow-scale resistance geometry is exactly symmetric under S3 at the
quotient level even though it emerged from the D12-labelled microscopic set.

## 4. Multiscale interpretation

The effective geometry has two nested structures:

    fast scale:
      three copies of C4;

    slow scale:
      one triangle connecting the three C4 sectors.

This is a concrete finite multiscale geometry produced by the energy landscape.

It is different from:
- the original Z12 cycle metric;
- the barrier ultrametric;
- the strict-operator resistance metric of ST230.

These objects should remain distinguished.

## 5. Gauge and scale

A global clock rescaling

    k_d -> rho k_d

rescales every resistance by

    1/rho

but leaves:
- the C4/triangle topology;
- all resistance ratios;
- the separation hierarchy;
- the normalized Laplacian eigenvectors

unchanged.

Thus the SHAPE of the two-scale resistance geometry survives the unsourced
global clock rate.

## 6. Refinement boundary

This two-level C4/triangle structure is not yet a refinement tower.

To become one, FIN would need a rule taking each effective node/cluster to a
new copy of the same or predictably transformed incidence structure.

No such self-similar recursion has been derived here.
