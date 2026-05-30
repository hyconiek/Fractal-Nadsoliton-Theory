# Strict alpha bounded shell-linear chi_11 source no-go

Status: `bounded-shell-linear-histogram-functionals-export-chi11-only-by-importing-d1-vs-d5-label`

## Linear functional summary

- Total weights: `729`
- Full-Aut invariant weights: `243`
- Pair-distinguishing weights: `486`
- Full-Aut pair-distinguishing weights: `0`
- chi_11-covariant weights: `54`
- Full-Aut invariant chi_11-covariant weights: `0`
- Allowed strict chi_11 source weights: `0`
- chi_11-covariant shell-label imports: `54`
- chi_11 support-size histogram: `{'2': 6, '3': 16, '4': 12, '5': 12, '6': 8}`

## Symbolic certificate

- `score_A1`: For A1/A11, h=[4,3,2,1,0,0], so L(A1)=4*w1+3*w2+2*w3+w4.
- `score_A5`: For A5/A7, h=[0,3,2,1,4,0], so L(A5)=3*w2+2*w3+w4+4*w5.
- `pair_difference`: L(A5)-L(A1)=4*(w5-w1); therefore any full-Aut invariant linear functional with w1=w5 cannot distinguish A1 from A5.
- `chi11_condition`: To transform as chi_11 on the two branch pairs, L(A5)=-L(A1) and L(A5) must be nonzero.
- `full_aut_no_go`: Combining w1=w5 with the pair-difference formula gives L(A5)=L(A1), incompatible with nonzero chi_11 covariance.

## Proof certificate

- `finite_domain`: All 3^6=729 weights w_i in {-1,0,1} for shell-linear histogram functionals are enumerated exactly.
- `symbolic_no_go`: L(A5)-L(A1)=4*(w5-w1); therefore any full-Aut invariant linear functional with w1=w5 cannot distinguish A1 from A5. Combining w1=w5 with the pair-difference formula gives L(A5)=L(A1), incompatible with nonzero chi_11 covariance.
- `enumerated_no_go`: Full-Aut invariant pair-distinguishing weights: 0; allowed strict chi_11 source weights: 0.
- `import_boundary`: There are 54 chi_11-covariant weights, and 54 of them have w1!=w5, i.e. they import the d1-vs-d5 shell label.
- `scope_limit`: This is exhaustive only for bounded shell-linear histogram functionals, not for all possible strict nadsoliton source mechanisms.

## Minimal chi_11-covariant examples

- weights=`[-1, 0, 0, 0, 1, 0]`, score_A1=`-4`, score_A5=`4`, imports_label=`True`
- weights=`[-1, 0, 1, 0, 0, 0]`, score_A1=`-2`, score_A5=`2`, imports_label=`True`
- weights=`[0, 0, -1, 0, 1, 0]`, score_A1=`-2`, score_A5=`2`, imports_label=`True`
- weights=`[0, 0, 1, 0, -1, 0]`, score_A1=`2`, score_A5=`-2`, imports_label=`True`
- weights=`[1, 0, -1, 0, 0, 0]`, score_A1=`2`, score_A5=`-2`, imports_label=`True`
- weights=`[1, 0, 0, 0, -1, 0]`, score_A1=`4`, score_A5=`-4`, imports_label=`True`

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.
- The result is a finite bounded shell-linear histogram no-go, not an exhaustive strict-source theorem.
- No QW-2191 discharge.
- No ToE closure.
