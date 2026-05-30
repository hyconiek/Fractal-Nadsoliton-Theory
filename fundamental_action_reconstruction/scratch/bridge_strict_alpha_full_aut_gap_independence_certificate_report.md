# Strict alpha full-Aut gap independence certificate

Status: `full-aut-d1-d5-d6-closure-has-independence-number-3-so-five-active-exact-cover-is-unsat`

## Finite model

- Ring: `Z_12`
- Full-Aut closed forbidden shells: `[1, 5, 6]`
- Target active count: `5`
- Maximum independent size: `3`
- Target five-active UNSAT: `True`

## Independence profile

- k=`0`: count=`1`, gap necklaces=`[]`
- k=`1`: count=`12`, gap necklaces=`[[12]]`
- k=`2`: count=`36`, gap necklaces=`[[2, 10], [3, 9], [4, 8]]`
- k=`3`: count=`16`, gap necklaces=`[[2, 2, 8], [4, 4, 4]]`
- k=`4`: count=`0`, gap necklaces=`[]`
- k=`5`: count=`0`, gap necklaces=`[]`

## Gap elimination certificate

- k=4 necklaces eliminated: `True`
- k=5 necklaces eliminated: `True`
- k=4 necklace count after d1 gap reduction: `8`
- k=5 necklace count after d1 gap reduction: `3`

## Proof certificate

- `graph_definition`: G has vertex set Z_12 and forbidden edges at folded distances d1,d5,d6, the full-Aut closure of the earlier d1+d6 exact-cover clauses.
- `gap_reduction`: For any independent support, the d1 ban forces every cyclic gap to be at least 2; all possible k-support gap certificates are therefore integer necklaces summing to 12.
- `four_support_elimination`: All 8 min-gap necklaces for k=4 contain a consecutive interval of folded length 5 or 6, so no 4-support survives.
- `five_support_elimination`: All 3 min-gap necklaces for k=5 contain a consecutive interval of folded length 5 or 6, so no 5-support survives.
- `tightness`: There are 16 independent 3-supports, while k=4 count is 0 and k=5 count is 0; hence alpha(G)=3.
- `selector_consequence`: The target exact-cover cardinality 5 is impossible after full-Aut clause closure; the successful d1+d6 selector remains chi_11-conditional.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, or unit-axis bit from strict nadsoliton geometry.
- The result is a finite full-Aut graph/gap no-go certificate, not a selector-origin theorem.
- No QW-2191 discharge.
- No ToE closure.
