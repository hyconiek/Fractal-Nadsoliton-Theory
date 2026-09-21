# R7N-044/R7N-046 full fixed-fixture phase theorem

## Exact fixture

The amplitudes are the exact declared decimal fixture

- r3 = 0.1131879146
- r4 = 0.1698528641
- r5 = 0.2269339093
- z6 = -0.3380663037

with normalized phase coordinates z = phi/(2*pi) on the full torus [0,1]^3 with periodic seams.

## Theorem

For this fixed full log-mgf phase fixture there are **exactly 60 critical points** on the full phase three-torus. Their negative-index histogram is

- index 0: 12,
- index 1: 24,
- index 2: 18,
- index 3: 6.

The 60 local roots have certified uniqueness collars. The smallest certified full collar radius is 3/10000 rad and the largest is 3/2000 rad. The minimum torus L-infinity separation of distinct root centers is 0.713167644237306 rad; the minimum separation after subtracting both collar radii is 0.712367644237307 rad, so all collars are pairwise disjoint.

## Complement proof

The torus complement is exhausted by a layered proof:

1. K16 baseline: 54,341 formula-certified gradient-exclusion leaves and 5,382 residual leaves.
2. K20 direct reclassification: 2,887 further gradient-exclusion leaves and 2,495 residual parents.
3. K20 adaptive closure: 22,405 gradient-exclusion leaves plus 864 leaves contained in certified root collars, with zero unresolved leaves.

Formula-level independent replay recomputed all 54,341 K16 leaves and all 25,292 K20 leaves: **79,633/79,633 passed, zero failures**. The independent geometry/mutation audit also passed.

## Relation to the quartic fixture

The separately certified quartic fixture also has exactly 60 critical points with the same index histogram. The catalogs carry a bijective label correspondence, but this campaign does **not** prove a global continuation/homotopy theorem between the two functions.

## Nonconclusions

This theorem is only for the stated fixed amplitudes and fixed negative alternating sign. It does not establish amplitude robustness, variable-amplitude phase exhaustion, a full-X7 stationary census, global minimizer uniqueness, or a physical coexistence/source law.
