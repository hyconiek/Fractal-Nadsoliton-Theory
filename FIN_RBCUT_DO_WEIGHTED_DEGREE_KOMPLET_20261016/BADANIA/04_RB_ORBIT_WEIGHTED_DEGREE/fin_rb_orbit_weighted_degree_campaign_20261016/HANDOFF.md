# FIN RB orbit flow: native ellipticity and weighted-degree criterion

Date 2026-10-16. Output `fin_rb_orbit_weighted_degree_campaign_20261016/`. Only new intrinsic mathematical research is packaged.

## Disposition

`PROVED_ACTUAL_MT_UNIFORM_NATIVE_EDGE_ELLIPTICITY_AND_EQUIVALENT_WEIGHTED_ORBIT_CUT_GATE; Q12_ALL_N_WEIGHTED_CUT_BOUND_OPEN`

**New actual FIN-specific theorem.** Because the real MT cylinder has all coordinates `p_j>103/1250` and A7 is the fixed PSD matrix with top eigenvalue Lambda6, every retained native opposite-parity count edge has `q0(x,y)>=a_0 N` with the **explicit conservative** `a_0=(103/15000)exp(-24)>0`, independent of N. The leave-one-out 12-way softmax and the exact rejection-reflected domain are retained.

Define `Q_N` as the largest positive exact drift-source sum on an orbit-union cut divided by `N` times its multiplicity-weighted number of outward **allowed** H edges. The full native source includes ALL reflected boundary states. Then for every finite N:

`Q_N <= C_N^vert <= K_N^* <= Q_N/a_0`, and `a_0 K_N^*<=C_N^vert<=K_N^*`.

Here `C_N^vert` is the smallest possible original exit-rate constant (total rate <= CN), and `K_N^*` is the smallest possible native reference edgewise cap on H. In particular **the actual original all-N vertex cap is equivalent to bounded `K_N^*` and equivalent to bounded `Q_N`**. This corrects the former *general-graph* warning in its application to the **local MT** reference; the old toy counterexample remains valid outside the uniform ellipticity assumption.

The exact weighted-degree score has an orbit representative formula needing at most 36 upward label-pair checks per orbit. It is not necessary to sum individual edge weights for this criterion. It is still necessary to establish an all-N orbit-cut inequality or exhibit real q12 countercuts; this executor did neither. No real q12 kinetic family, no time-integrated IBA obstruction, no physical activity source.

## Files

- `UNIFORM_NATIVE_EDGE_ELLIPTICITY.md` — model-specific occupation floor, softmax lower bound, explicit `exp(-24)` constant, Gibbs neighbour comparability.
- `WEIGHTED_ORBIT_CUT_EQUIVALENCE_THEOREM.md` — exact theorem, converse, orbit-representative formula.
- `ACTUAL_Q12_NEXT_GATE.md` — precisely unpaid all-N weighted-cut inequality / countercut.
- `NEGATIVE_CONTROLS.md` — necessity of locality and separation from generic examples.
- `verify_weighted_degree.py`, `FINITE_CHECKS.json` — exact rational finite controls.
- `REPORT.md`, `RESEARCH_SUMMARY_PL.md`, `CLAIM_REGISTER.json`, `RESULTS.json`, `SOURCE_PROVENANCE.json`, `REPLAY.md`, `INPUTS.sha256`, `MANIFEST.sha256`, `PROPOSED_AGENTS_PATCH.md`.

## Source limitations

The actual archived MT-001 source and current local preceding handoffs were examined. This environment did not expose a verified current remote Git checkout, the latest RB/FP review sources or a full original versioned A7 `safe_seed.py` as independent accessible inputs. The proof uses only the inherited *structural* PSD / top eigenvalue and LOO contracts, **not a guessed numerical A7 matrix**. Source hashes are bound separately; synthetic replay does not pretend to validate current custody. Previous files, plans, AGENTS.md and primary FIN are unchanged.

## Stop

Do not promote a bounded actual-q12 `Q_N` until proved, do not claim unconditional IBA failure or a unique physical kinetic source, and do not change any QM/GR/SM/ToE gate. Next valid task: prove or refute the single weighted orbit-cut inequality for genuine MT q12/A7.
