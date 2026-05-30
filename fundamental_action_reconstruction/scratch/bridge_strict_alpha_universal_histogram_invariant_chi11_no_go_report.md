# Strict alpha universal histogram-invariant chi_11 no-go

Status: `full-aut-invariant-histogram-only-sources-cannot-distinguish-a1-from-a5-or-export-chi11`

## Histogram orbit summary

- Histogram classes: `35`
- Full-Aut histogram orbits: `24`
- Fixed histogram orbits: `13`
- Two-histogram orbits: `11`
- Invariant Boolean classifier total: `16777216`
- Singleton A5-not-A1 invariant classifiers: `0`

## A1/A5 orbit certificate

- A1 histogram: `[4, 3, 2, 1, 0, 0]`
- A5 histogram: `[0, 3, 2, 1, 4, 0]`
- swap_d1_d5(A1): `[0, 3, 2, 1, 4, 0]`
- Same full-Aut histogram orbit: `True`
- Orbit support count: `24`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are enumerated and collapsed to distance-histogram classes.
- `histogram_swap_fact`: A1 has histogram [4,3,2,1,0,0] and A5 has [0,3,2,1,4,0], exactly the d1<->d5 swap.
- `universal_invariant_no_go`: For any histogram-only source F invariant under full Aut, F(h)=F(swap_d1_d5(h)); hence F(A1)=F(A5) and F cannot export chi_11.
- `classifier_count`: The 35 histogram classes form 24 full-Aut histogram orbits; invariant Boolean classifiers are constant on these orbits, so singleton A5-vs-A1 classification count is 0.
- `anti_invariant_boundary`: A histogram expression can flip between A1 and A5 only if it is anti-invariant under d1<->d5, which imports the missing shell-label orientation.
- `scope_limit`: This no-go is universal for histogram-only full-Aut invariant sources, not for all possible non-histogram strict nadsoliton sources.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.
- The result is a universal histogram-only full-Aut invariant no-go, not an exhaustive strict-source theorem.
- No QW-2191 discharge.
- No ToE closure.
