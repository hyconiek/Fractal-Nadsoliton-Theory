# Strict alpha quadratic histogram invariant chi_11 no-go

Status: `degree2-full-aut-histogram-invariants-cannot-distinguish-a1-from-a5-or-export-chi11`

## Summary

- Invariant orbit sums: `21`
- Invariant pair-distinguishing orbit sums: `0`
- Anti-invariant orbit differences: `6`
- Anti-invariant chi_11-covariant differences: `5`
- Allowed strict chi_11 source count: `0`

## Symbolic certificate

- `histogram_relation`: A1/A11 has h=(4,3,2,1,0,0) and A5/A7 has h'=(0,3,2,1,4,0), exactly the d1<->d5 swap of h.
- `invariant_principle`: Any full-Aut invariant histogram polynomial P satisfies P(h)=P(swap_d1_d5(h)), so P(A1)=P(A5).
- `degree2_enumeration`: The degree<=2 monomial basis has 6 linear + 21 quadratic monomials; quotienting by d1<->d5 gives 21 invariant orbit sums.
- `anti_invariant_boundary`: Orbit differences can flip sign under d1<->d5, but choosing an orbit difference chooses an orientation of the d1/d5 pair and therefore imports the missing shell-label bit.

## Proof certificate

- `finite_domain`: All degree<=2 histogram monomials and their d1<->d5 full-Aut orbits are enumerated exactly.
- `invariant_no_go`: Every one of the 21 invariant orbit sums has equal value on A1/A11 and A5/A7, so their linear span cannot export chi_11.
- `symbolic_reason`: Any full-Aut invariant histogram polynomial P satisfies P(h)=P(swap_d1_d5(h)), so P(A1)=P(A5).
- `anti_invariant_boundary`: Five anti-invariant orbit differences are nonzero chi_11-covariant, but each imports an orientation of the d1/d5 shell pair.
- `scope_limit`: This is exhaustive only for degree<=2 histogram invariants, not for all possible strict nadsoliton source mechanisms.

## Anti-invariant chi_11-covariant rows

- `d5 - d1`: A1=`-4`, A5=`4`, imports_label=`True`
- `d5*d5 - d1*d1`: A1=`-16`, A5=`16`, imports_label=`True`
- `d2*d5 - d1*d2`: A1=`-12`, A5=`12`, imports_label=`True`
- `d3*d5 - d1*d3`: A1=`-8`, A5=`8`, imports_label=`True`
- `d4*d5 - d1*d4`: A1=`-4`, A5=`4`, imports_label=`True`

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.
- The result is a finite degree<=2 histogram invariant no-go, not an exhaustive strict-source theorem.
- No QW-2191 discharge.
- No ToE closure.
