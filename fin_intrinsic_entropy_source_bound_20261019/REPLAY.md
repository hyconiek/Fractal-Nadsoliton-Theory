# ID-001 — reproduction instructions and scope

1. Examine `SOURCE_CLASS.json`, the two analytical proofs, and the controls.
2. From this handoff folder run `python3 verify_id001.py --output results_replay.json` (requires standard Python, NumPy; no network or original FIN scripts). Compare relevant numbers with `results.json` using floating tolerance. This diagnostic reconstructs A7 from the displayed k(d), verifies the rank-7 Fourier row by a separate matrix multiplication, computes p at **one fixed g=3**, checks row sums, stationarity, two entropy formulas and independent mean jump-rate equality. It neither checks KS entropy of a flow nor proves the entropy inequality.
3. Check algebra by hand: `hdelta = -sum p_i [(1−t(1−p_i))log(1−t(1−p_i)) + sum_{j!=i} t p_j log(t p_j)]` with `t=1−exp(−2delta)`. The coefficient of `delta log(1/delta)` is `rbar=2(1−sum p²)`, not 2 or `2sum p_i`.
4. Independently replay the block-chain-rule inequality from `ENTROPY_BUDGET_PROOF.md` and test the monotonicity range of Fano inverse (`e≤11/12`).
5. Verify source custody only for bytes actually copied into this environment; `INPUTS.sha256` is not a root Git manifest. `python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only` cannot succeed on the partial extracted archive: absent root scientific inputs are a source-availability limitation. Run it in the user's actual repository after checking that HEAD and scientific inputs have not changed.
6. Compare generated output against `MANIFEST.sha256` using `sha256sum -c MANIFEST.sha256` from inside this result folder. Manifest omits itself and temporary `results_replay.json` by design.

No heavy computation, long parameter sweep, stochastic trajectory simulation, source mutation, or external communication was carried out.
