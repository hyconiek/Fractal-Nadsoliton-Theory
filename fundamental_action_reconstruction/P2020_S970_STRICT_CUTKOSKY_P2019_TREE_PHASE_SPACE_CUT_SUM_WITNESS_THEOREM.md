# P2020 S970 Strict Cutkosky P2019 Tree Phase-Space Cut-Sum Witness Theorem

Status: `P2020_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_EXPORTED_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

P2019 exported the first strict-side local transverse tree component for

```text
graviton -> gauge_gauge.
```

P2020 performs the next conservative step: integrate that component over the
exact massless two-body phase space.  In response to the main rigor risk in the
previous version, P2020 now first rederives the P2019 axis-frame amplitude in a
generic spherical real transverse-polarization frame instead of merely assuming
angle-independent transport.  It does not add a new amplitude ansatz and it does not promote the
result to dressed Cutkosky closure.

## Input from P2019

P2019 gives the local tree-level transverse amplitude matrices, in units of
`kappa*Z_gauge`,

```text
M_plus  = [[-2, 0], [0, 2]],
M_cross = [[ 0,-2], [-2,0]].
```

Therefore the local polarization sum is

```text
sum_pol |M_tree_transverse|^2 / (kappa^2 Z_gauge^2) = 16.
```

## Generic angular-transport certificate

Let

```text
k1 = (1, n),
k2 = (1,-n),
n = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)).
```

P2020 constructs the real transverse-polarization basis

```text
e_theta = (0, cos(theta)cos(phi), cos(theta)sin(phi), -sin(theta)),
e_phi   = (0,-sin(phi), cos(phi), 0).
```

The script verifies symbolically that these vectors are transverse and
orthonormal for both `k1` and `k2`, then rebuilds the same Maxwell cross-stress
contraction used in P2019.  The generic-frame matrices are again

```text
M_plus  = [[-2, 0], [0, 2]],
M_cross = [[ 0,-2], [-2,0]],
```

so the polarization sum `16` is no longer only an axis-frame value; it is a
symbolically transported real transverse-polarization-frame result for this tree
component.  P2020 also keeps the external graviton linear-polarization indices before taking any trace:

```text
G_{P P'} = sum_{a,b in {theta,phi}}
  M_{P,ab} M_{P',ab}

G = [[8,0],[0,8]]  in the real {plus,cross} linear-polarization basis.
```

This is intentionally **not** a circular helicity-eigenstate claim; the previous
wording was too strong.  P2020 now describes the object as
linear-polarization-resolved.

## Exact phase-space measure

P2020 then uses the standard massless two-body center-of-mass measure

```text
dPhi2 = (1/(32*pi^2)) dphi dx,
phi in [0,2*pi], x=cos(theta) in [-1,1].
```

Because the generic angular-transport certificate rederives an angle-independent
polarization sum,

```text
CutSum_tree_linear_polarization_matrix_no_identical_symmetry
  = integral dPhi2 * G
  = [[1/pi,0],[0,1/pi]]

trace(CutSum_tree_linear_polarization_matrix_no_identical_symmetry) = 2/pi
```

and, if the two gauge bosons are treated with an identical-final-state symmetry
factor `1/2`,

```text
CutSum_tree_linear_polarization_matrix_identical_final_state
  = [[1/(2*pi),0],[0,1/(2*pi)]],

trace(CutSum_tree_linear_polarization_matrix_identical_final_state) = 1/pi.
```

Thus P2020 exports both the linear-polarization-resolved matrices and the scalar trace
convention window

```text
[1/pi, 2/pi]
```

in units of `kappa^2 Z_gauge^2`.  The scalar entries are traces of the matrix
object, not a replacement for the linear-polarization-resolved cut-sum.

## Numerical cross-check

The script also runs `scipy.integrate.dblquad` over the same domain and compares
against the exact SymPy result.  The quadrature is a verification cross-check;
the theorem object is the exact symbolic integral.

## No-false-pass boundary

P2020 is still only a tree-component cut-sum witness.  The angular-transport
certificate removes one local-frame weakness, but it does **not** supply loop
dressing, BRST cohomology, or a discontinuity calculation.  It is not:

1. a dressed loop amplitude,
2. a BRST cohomology theorem,
3. a same-scheme `MSbar_B1_seed` dressed object,
4. a `DiscM = CutSum` proof,
5. all-state unitarity,
6. `QW-2191` discharge,
7. ToE closure.

The P1953 contract is therefore updated only partially; even the new
linear-polarization-resolved matrix remains tree-level and not dressed:

```text
phase_space_measure and integration_domain = PARTIAL_TREE_COMPONENT_EXACT_DPHI2_EXPORTED
AbsM_dressed_squared_common_basis = PARTIAL_TREE_TRANSVERSE_POLARIZATION_SUM_EXACT_NOT_DRESSED
CutSum_common_basis = PARTIAL_TREE_LINEAR_POLARIZATION_RESOLVED_CUTSUM_CONVENTION_WINDOW_NOT_DRESSED_NOT_DISC_MATCHED
```

## Next honest step

Derive the matching same-scheme `DiscM_common_basis` or a first dressed-residue
factor for exactly this normalized linear-polarization-resolved tree cut-sum component.  Only after both
sides are in the same scheme and normalization should the project test
`DiscM_minus_CutSum_simplified`.
