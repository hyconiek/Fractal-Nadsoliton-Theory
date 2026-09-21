# FIN rank-seven — continuation checkpoint 2026-09-15

**Start with `HANDOFF.md`.**

This directory is the 2026-09-15 continuation of the completed 2026-09-14 128-task campaign. The original handoff is preserved as `HANDOFF_BASELINE_20260914.md`. The R7P ledger remains frozen; new research is tracked separately in `CONTINUATION_LEDGER.json` as FR frontier work.

The current scientific target is the still-open residual physical positive-orthant four-amplitude ceiling. The continuation adds global projected tails and a family of interval-certified anisotropic local boxes, but it does **not** claim a complete global 4D theorem.

---

# fin_rank7_followup

Execution package for `FIN_Post_Handoff_Research_Master_Plan_EN.md`. Imported inputs under `inputs/` are immutable research evidence.

## Current scientific state

Work package G (R7P-041--056) is closed on its exact boundary-Ising domain:
`lambda2(Cov)<=sigma_*` throughout the complete four-state boundary closure,
with unique double-root equality.

Work package H (R7P-057--064) is now also closed on its stated shared-field
domain.  Exact/interval-certified complementarity gives

`lambda2(W_par)<=sigma_*` for all `J3,J4,J5,J6>=0`.

The conservative proof uses a freshly certified `p_dom>711/1000`,
`lambda2(C_+)<511/2000`, and a positive final Weyl margin.  The exact numerical
dominant-mass minimizer is retained only as a numerical reproduction.

This is **not** yet the full off-face four-amplitude theorem for
`M=W_par+b b^T`, and is not a full-seven-coordinate theorem.  Work package I
(R7P-065--072) remains the next core lane.

## Packaging note

The supplied outer handoff archive is missing
`inputs/FIN_research_artifacts_pre_and_post_Discord.zip`, although its hash
remains frozen in `INPUT_HASHES.json`. Consequently top-level `verify.py`
correctly reports this historical packaging defect. The expected hash has not
been changed to hide it.

## Work package I continuation — 2026-09-14

R7P-065 is closed: the full four-amplitude target is frozen using a stable
Schur/inertia representation with a globally positive scalar denominator
`eta=1-lambda6 q(1-q)/(3 sigma_*)`.

R7P-066 is a numerical-only adversarial search ledger. Three fixed-seed DE
runs over the **entire exact compactified closure** plus structured slices find
no positive `lambda2-sigma_*` beyond roundoff at the known equality point.
This is not a proof.

R7P-067 is exact: the whole nonnegative field orthant has a compact closure.
In the chart
`r=e^-2J3, s=e^-3J4/2, t=e^-J5/2, y=e^-2J6`, the `y=0` face is exactly the
already-certified G boundary cube, while odd weights retain the common
`t^(2+-sqrt(3))` dependence. No finite-field tail cutoff is needed.

R7P-068 is in progress. The first-order physical cone is safe and the formerly
numerical reoptimized parity slope is now interval-certified near
`-0.1312828584000`, but an explicit finite 4D radius still requires a validated
mixed higher-order remainder.

Current tests: `68 passed, 1 deselected` when excluding the one mutation test
whose baseline necessarily fails because the supplied package lacks the frozen
historical input ZIP. Full `verify.py` still reports only that pre-existing
missing-input defect.


## Final campaign handoff (2026-09-14)
The 128-task campaign has terminal dispositions for all tasks. Start with `HANDOFF.md`, `REPORT_FINAL.md`, `verification_final.json`, and `REPLAY.md`. Scientific completion is scoped: unresolved research frontiers are preserved explicitly.

## Post-FR17 checkpoint

The latest portable checkpoint includes durable FR5, FR15 and FR16 certificates plus the FR17 fixed-seed compact-residual search. `POST_FR17_MANIFEST.sha256` is the authoritative hash ledger for this updated checkpoint; historical manifests are retained for provenance.
