# FIN RG-E01–RG-E05 — full bounded operational identifiability handoff

**Date:** 14 October 2026. **Architect branch:** `SOURCE_AND_OPERATIONAL_IDENTIFIABILITY_GATE`. **Executed:** RG-E01–RG-E05 only. **RG-E06:** not executed (architect-only).

## Executive verdict

`DESIGN_ONLY / NO_NEW_SOURCE / NO_NEW_PHYSICAL_INPUT`, with mathematical model theorems and independently reproduced finite numerical design contrasts. The full parent is a mathematical reference; no real channel measurement or empirical validation was performed. All `model/physical-premise/bridge/prediction/evidence` roles are recorded.

## 1. Actual input contract

Current user-attached operational plan and the complete predecessor RG-D handoff ZIP are preserved verbatim under source_context, with SHA-256 in the manifest. Current repository files were reviewed via public raw pages (`PHYS-002`, `PHYS-004`, `PHYS-005`, `PHYS-014`, `RG_D06_DECISION`, `AGENTS.md`), but a full byte-verifiable git checkout/HEAD was not available; source code is a declared independent reimplementation, not a reproduced repository executable. Prior user-pasted intake is treated as user-supplied scientific context, not an independently audited original byte file.

## 2. Main new scientific results

| Fixed ideal hypotheses | N2 difference TV at g=3 | M3 ordered TV at old g=5.145... |
|---|---:|---:|
| FIN vs full Potts | 0.096998787035 | 0.085011536357 |
| FIN vs flat P7 | 0.011787191579 | 0.010157172973 |
| Full Potts vs flat P7 | 0.106340140736 | 0.094052748124 |

N2 *static* TV bounds were recomputed from exact-decimal strict-kernel inputs with an independent 65-decimal-digit interval implementation. The 1e-15-scale **public outward endpoints are conditional on the mpmath.iv backend** and a supported fixed finite model, not full IEEE proof. The M3 numbers are numerical comparisons only; the two parents must not be conflated.

Full M3 conditional log-odds contrasts exactly identify only a centered scaled interaction quotient when hidden b/c conditioning exists; this observation is NOT available in the PHYS-014 N2 proposal. Clock and gain-scale freedoms remain. FIN nonzero Fourier ratios relative to mode k3: 1.000000000, 1.121424061, 1.171917116, 1.194133704. Flat P7 ratios are all one and full-Potts has extra positive low modes, under the specified fixed models.

For completely uncalibrated detector readout, a constant-output stochastic channel proves minimax separation zero. Under *hypothetical* ±1% gain and ≤0.005 detector-TV design allowances, the conservative FIN/flat-P7 guarantee is zero. With no gain or detector uncertainty, the ideal FIN/flat P7 difference is 0.011787191579; half of it (0.005893596) is the equal per-model TV-error separation budget. Not a claim that the hardware achieves it.

N2 kinetics with **same supplied equilibrium law** and model-g0-normalized rates yield distinct *two-time* 13×13 joint laws at model t=0.4 and t=1.0, starting from g0 uniform preparation and switching to g=3. At FIN, TV(heat-bath, Metropolis)=0.256627921020, TV(heat-bath, Barker)=0.194876586235. The ideal invalid outcome is retained explicitly with zero mathematical probability; its **physical** rate is unknown.

## 3. Verification / deliverable scope

- Independent N2 generator built by a different dense implementation and FIN X X^T rather than the producer's circulant shortcut: **58/58 controls pass**. The reference propagation uses dense `scipy.linalg.expm`; producer uses sparse `expm_multiply`.
- Separate N2 interval producer (mpmath.iv conditional) outputs probability and pairwise TV endpoint enclosures; comparison claims are exact-model, not experiment-level.
- Each task RG-E01–RG-E05 contains `REPORT.md`, `RESULTS.json`, `NONCONCLUSIONS.md`, `REPLAY.md`, `INPUTS.sha256` and its required special artifacts. Script `run_operational.py`, source `source_model.py`, static independent interval `interval_static.py`, `independent_verify.py` and complete output arrays are included. Prior RG-D is embedded unmodified.
- The local original source byte/input hashes are recorded; public source pages viewed via web are **not** represented as locally verified original SHA.

## 4. Final boundary and one executable next step

`RG-E05/PREREGISTRATION.json` is a prospective **N2 engineered** design only, not an approved experiment. It requires independent PHYS-014 actual 24-channel calibration/weights, physical g0 clock, detector/error/invalid logs and disjoint validation before empirical status can be discussed. The separate independent A7 source atom is still missing. Architect RG-E06 must decide whether to retain `PHYSICAL_TEST_DESIGN_ONLY` / `NO_NEW_PHYSICAL_INPUT` or authorize a separately scoped operational measurement; executor has made no architectural decision.

Do not open task335, redo B4/FIN-M0, invent new dimensions/particles, or infer QM/GR/ToE.

See `CLAIM_REGISTER.json`, `NONCONCLUSIONS.md`, `REPLAY.md`, `SOURCE_ACCESS.md`, `RG-E03/FALSIFICATION_TABLE.md`, `RG-E04/ERROR_BUDGET.md`, and `RG-E05/STOP_CERTIFICATE.json`.
