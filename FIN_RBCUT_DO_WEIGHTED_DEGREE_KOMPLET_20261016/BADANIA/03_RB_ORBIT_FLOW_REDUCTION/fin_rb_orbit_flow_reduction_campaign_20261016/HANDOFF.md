# FIN RB-CUT: Orbit flow reduction handoff

2026-10-16. New scoped research package. **Disposition: `EXACT_FINITE_ORBIT_REDUCTION_ACCEPTABLE_FOR_REVIEW; ACTUAL_ALL_N_Q12_CUT_AND_VERTEX_CAP_OPEN`.**

## What is new

1. A rigorous specialization of finite group averaging to the *actual* q12 parity-preserving dihedral subgroup. Full finite original graph feasibility for native-relative edge caps is equivalent to feasibility on the exact orbit network, with a constructive lift. Consequently sharp K_N^* is decided by orbit-union cuts.
2. A second equivalence, not confined to dominated rates: the actual total count exit bound `sum_y q(x,y)<=CN` is also exactly preserved by the orbit flow quotient. This reaches the original uniform vertex-cap class.
3. A complete finite Farkas obstruction (`phi,beta`) for this quotient vertex-cap LP. It distinguishes capacity-free lower-ideal obstructions from capacity bottlenecks.
4. Separate rational tests and dual/primal certificates, including a non-feasible graph control and complete full-cut enumeration for a six-vertex symmetric toy.

## What is NOT proved

No all-N true q12 inequality, no valid actual q12 negative orbit ideal, no uniformly capped actual diagnostic generator, no temporal IBA failure or universal IBA theorem. Group reduction reduces size but does not provide analytic control of irregular MT shell cuts. Any result of a toy fixture is explicitly non-FIN.

## Files

- `ORBIT_QUOTIENT_AND_LIFT_THEOREM.md` — precise hypotheses, finite proofs, explicit lift, dual.
- `ACTUAL_Q12_NEXT_GATE.md` — fully typed remaining actual all-N problem.
- `EXACT_RATIONAL_FIXTURES.md`, `verify_orbit.py`, `EXACT_CHECKS.json` — independently replayable exact rational controls.
- `REPORT.md`, `RESEARCH_SUMMARY_PL.md`, `CLAIM_REGISTER.json`, `NONCONCLUSIONS.md`, `REPLAY.md`, `INPUTS.sha256`, `MANIFEST.sha256`.
- `PROPOSED_AGENTS_PATCH.md` — only a proposal; no modifications of AGENTS.

## Handoff instruction

Attempt just one analytically grounded estimate for the **true** quotient orbit network, with the original MT domain, departure set, fixed phase threshold, exact Gibbs weights and reflected HB source. A good outcome is an all-N bound, a true dual obstruction with correct N quantifiers, or the exact remaining inequality. Do not replace the actual q12 graph with the toy or infer that failure of edgewise native domination rules out the original vertex-cap class.
