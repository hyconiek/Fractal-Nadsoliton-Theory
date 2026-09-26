# REVERSIBLE-CONDUCTANCE-REFINEMENT-BRIDGE-115
## The Dirichlet conductance of the reversible heat-bath network obeys exactly the FIN series-refinement algebra

Date: 2026-09-26

Status:
- exact reversible-network theorem;
- mathematical typed bridge between metastable dynamics and the static
  resistance algebra;
- NOT an identification with physical stiffness or physical distance.

## 1. Reversible metastable generator

For a reversible effective Markov chain with stationary weights pi_i and rates

    k_ij,

define the symmetric edge conductance

    C_ij = pi_i k_ij = pi_j k_ji.

Its Dirichlet form is

    E(f)
      = (1/2) sum_(i,j)
          C_ij (f_i-f_j)^2.

This is the natural quadratic energy associated with the reversible generator.

## 2. Exact three-node series elimination

Take a chain

    A -- B -- C

with conductances

    C1, C2

and no direct A-C edge.

For fixed endpoint values f_A,f_C, minimize

    E
      = C1(f_A-f_B)^2
        +C2(f_B-f_C)^2

over the internal value f_B.

The minimizer is

    f_B
      = [C1 f_A+C2 f_C]/[C1+C2].

Substitution gives

    E_eff
      =
      C_eff (f_A-f_C)^2,

where

    boxed:
    C_eff
      =
      C1 C2/(C1+C2).

Equivalently,

    boxed:
    1/C_eff
      =
      1/C1+1/C2.

This is exactly the same algebra as the accepted FIN static edge composition

    c(a+b)
      =
      c(a)c(b)/[c(a)+c(b)].

## 3. General statement

For a reversible network, eliminating internal nodes from the Dirichlet form is
Kron/Schur reduction of the weighted graph Laplacian.

Series edges therefore compose by additive resistance.

So the relevant bridge object is NOT barrier height.

It is the reversible Dirichlet conductance.

## 4. Relation to the maximum-entropy heat bath

Reports 54-57 derived an exact reversible heat-bath/Gibbs process and its
Onsager/Dirichlet structure.

Reports 110-112 then reduced its localized regime to an effective reversible
twelve-minimum chain.

Therefore the metastable chain supplies, without inventing a new algebra, a
conductance matrix

    C_meta.

Its Schur elimination obeys the same resistance composition law already used
in FIN refinement.

This is a genuine algebraic compatibility between two previously separate
research lanes.

## 5. Clock gauge

Multiplying the microscopic refresh rate by rho gives

    k_ij -> rho k_ij,

hence

    C_ij -> rho C_ij,

and

    R_ij=1/C_ij -> R_ij/rho.

Thus all resistance RATIOS and normalized resistance geometry are invariant
under the global clock gauge.

Only one overall scale remains free.

This closely parallels the static refinement result

    c(ell)=kappa0/ell,

which also leaves one global positive scale kappa0.

## 6. What is still missing

The compatibility does NOT prove

    metastable resistance = physical length.

Missing typed arrows include:

1. why the effective twelve-minimum chain is a physical incidence graph;
2. why its Dirichlet conductance is the same physical object as the static
   edge coefficient c;
3. how one metastable cell refines into multiple cells;
4. how the state-dependent rates depend on a geometric length ell;
5. any SI length or clock calibration.

So the result is:

    SAME COMPOSITION ALGEBRA,

not:

    SAME PHYSICAL QUANTITY.
