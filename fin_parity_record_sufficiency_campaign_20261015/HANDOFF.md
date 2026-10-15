# FIN LI-001–LI-005: bounded full-parity predictive risk — handoff

**Date:** 2026-10-15. **Verdict:** `CONDITIONAL_BOUNDED_PARITY_RISK_THEOREM_SUBMITTED_FOR_INDEPENDENT_INTAKE`.

**Scope:** new LI work only. Earlier RI/UF/MT/OBS packages are *source references* bound in `SOURCE_INPUTS.sha256`, not included or repeated in this handoff. The original labelled FIN dynamics, A7, gain, unit clock, domain, and physical guardrails remain unchanged.

## Main result (conditional analytical theorem)

For `N∈12ℕ`, q12 labelled LOO-HB at exact `g_N=g_c`, iid-uniform t=0 and the deterministic post-MT local window, let `O_N` be the complete **one-tag parity** trajectory through scaled time 1/2, and `C_N=sign(Z_N(1))` with the three finite-N outcomes `−,0,+`. Under the inherited MT local-window and UF/OBS stationary critical-path assumptions, plus the new LI likelihood and conditional-cylinder proof,

`r_N(O_N→C_N) → r_*(Y_[0,1/2]→sign Z(1))`,

where `dZ=−Z³/3 ds+√2 dB`, `dY=Z ds+√2 dW`, `Z(0)~ρ_0∝exp(−z⁴/12)`, and B,W are independent. Risks use 0–1 loss, not entropy or total variation between all observation laws.

### New LI work

- **LI-001:** independently verified exact global path-TV feedback bound `TV(P,R)≤|g|M_A√((N−1)T)/(2N)` for equal initial laws and every T, and derived a separate local-Gibbs initial-parity/product and target projection payment. This bound covers model time through the *future target* `T=√N`, not merely the observation time.
- **LI-002:** exact sign-flip (not attempt-mark) likelihood `log L_N=∫log(1−S_-q_N)dJ+(1/2)∫Sq_Ndt`, exact score `(1/2)∫q_NdM_S`, and time-integrated high-probability error using MT stationary-reflected moments plus path-event TV transfers.
- **LI-003:** joint fair reference background/record/score/future-sign limit. Conditional Gaussian unit-mass identity, together with finite likelihoods' exact mean one, pays **uniform integrability**, not only weak convergence.
- **LI-004:** full finite-cylinder conditional-density and bounded Bayes-risk limit. The crucial estimate approximates `1_{C=c}L_N` in L1 by bounded continuous functions of finitely many background/observation coordinates, integrates out background under the paid product Q reference, and then refines cylinders *after* N→∞.
- **LI-005:** strict source/status/physical-scope handoff; no automatic follow-on campaign.

### What must be independently audited

The submitted proof invokes accepted UF/MT/OBS/RI mathematical input contracts. The mounted prior ZIPs contain authorial packages, but the named 2026-10-15 review/proof-completion files (L02/L03/L04 and certain MT/OBS reviews) were **not** available as actual review artifacts. These unavailable bytes are not hashed or claimed checked. The load-bearing checks are: (a) full labelled post-MT Gibbs/reflect coupling and boundary payment; (b) transfer of stationary-reflected `∫w²` and `∫m⁶` path-events through P→R and the product reference without multiplying an o(1) TV error by √N; (c) joint stopped stochastic integral convergence; (d) finite-cylinder L1 conditional projection.

If those inherited items fail independent intake, retract the dependent end-to-end risk theorem, **not** the exact finite feedback KL/TV theorem or parity likelihood identities.

### Reproduction

`independent_checks.py` runs finite arithmetical source-convention checks against an unrelated small q=3 synthetic kernel; it does *not* certify q12 A7 geometry, process limits or the Bayes theorem. Run it before `verify_package.py`; inspect `CLAIM_REGISTER.json` and per-task proofs. Hash-binding via `SOURCE_INPUTS.sha256`, archive files via `MANIFEST.sha256`.

**Next architect checkpoint:** independent review of the existing source proof completions and this finite-cylinder/occupation derivation, not a new physical, entropy, computational-precision or kernel campaign.
