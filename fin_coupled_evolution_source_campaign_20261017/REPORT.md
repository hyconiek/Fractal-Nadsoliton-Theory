# FIN UE-001–UE-004 — bounded source-bound coupled-evolution campaign

Date: **2026-10-17**. Scope: analytical model-theory research and provenance audit, **not an independent scientific review or a physically validated unified law**.

## Terminal scientific verdict

**`CONDITIONAL_COUPLED_FIELD_ONLY_WITH_NEW_NONCOMMUTING_BOUNDS_TRANSFER_WITNESS_AND_REDUCTION_OBSTRUCTION; ACTUAL_SOURCE_BINDING_OPEN`**.

The requested next theoretical step was **not** narrowed to a classifier, q12 Hall transport, old far-source inventory, or a new phenomenological particle law. We examined the *declared existing action/operator interface* and extracted three precise reusable model facts, with transparent countercontrols. We cannot certify they hold for actual FIN coefficients because the 2026-10-17 audit source directory, seven audit files, QW scripts/reports and action contract were not retrievable in this run. No guessed array stands in for them.

## UE-001 — source/action

The reported 13-real-coordinate QW polynomial family is not the reported 12-complex-plus-1-real full action. The former has 26 canonical mechanical coordinates, the latter 50 real canonical field coordinates per spatial point before constraints. A real slice is invariant only if the full conjugation-symmetry and uniqueness hypotheses hold. The requested `M2`, `K_s`, `Y`, nonlinear coefficient arrays, signature and Cauchy domain have **not been bound**. The actual coefficient source remains unproved. `SOURCE_DELTA.md`, `PREMISE_LEDGER.md`, `ACTION_INTERFACE.json`, `SOURCE_ACCESS.md` record this explicitly.

## UE-002 — genuinely noncommuting form bounds

For the conditional scalar-kink family `Phi_K=v tanh(kx)`, with symmetric `A=M2+K_s`, `B=2v²Y`, the exact noncommuting identity is `V(x)=V_inf-B sech²(kx)` and `V_inf=A+B`.

- **Positive criterion:** `A>=0` and `V_inf>=0` imply `L_perp>=0` by convex interpolation.
- **Quantitative lower bound:** `inf spec L >= mu-k²s²`, where `mu=lambda_min(V_inf)`, `b=max(0,lambda_max(B))`, `s=(sqrt(1+4b/k²)-1)/2`.
- **Negative-direction criterion:** if `lambda_min(V_inf-2B/3)<-k²/3`, a trial `u sech(kx)` proves `inf spec L<0`, with no commutativity assumption.
- **Adverse exact control:** a real symmetric 2×2 **noncommuting** potential has far-vacuum eigenvalues `{1/2,3/2}` but Rayleigh quotient `-4/3`, certifying a discrete negative mode. **These matrices are NOT the actual FIN coefficient arrays.**

This answers a mathematically substantive part of transverse kernel-sensitive stability without pretending to compute a missing 24-real-component action Hessian. See `STABILITY_OR_NEGATIVE_DIRECTION.md`.

## UE-003 — exact conditional nonlinear transfer for nonzero Psi

For any actual canonical coupled action that admits `V=U+W+V_int`, the three sector-energy work identities follow directly from its own equations. Under local wellposedness, nonzero interaction gradients `F=d_Phi V_int`, `G=grad_Psi V_int`, and specially chosen **admissible initial velocities**, energy is transferred from the Phi diagnostic sector to Psi over a short positive interval while `dE_int/dt|_0=0`. This is exact under its hypotheses, with no new driving/bath/controller. A synthetic scalar case gives `E_Phi'(0)=-16delta/15`, `E_Psi'(0)=+16delta/15`. The same-field two-Gaussian profile has an exactly computed nonzero interaction-energy cross term `sqrt(2*pi) exp(-2a²)`; this is not a new two-body coupling. No generic scattering, spontaneous emission or two-wall absorption is claimed.

## UE-004 — quantified coarse closure obstruction

For the *frozen* projection `C(Phi,Phi_t,Psi,Psi_t)=(Phi,Phi_t)`, equal coarse initial states with different hidden accelerations (L² size `a>0`) yield separation at least `a t²/4` for an explicitly bounded small-time interval; any autonomous coarse prediction incurs worst-case error at least `a t²/8`. The same synthetic polynomial interaction gives `a=4epsilon²/sqrt(15)`. This is a rigorous conditional **non-closure result**, not proof of a strict fractal hierarchy. A memory-bearing reduction remains open, and no spatial refinement/composition map was supplied.

## Reproducibility and access limits

`python3 REPLAY/verify.py` executed **11 exact analytical algebra/integral controls** successfully. The code and the detailed results are included; no QW-program runs or full eigenspectrum campaign were performed. All unavailable source files are reported, not replaced with invented hashes. One materialized older AGENTS Library snapshot is hashed and retained for audit, clearly distinguished from current upstream.

**Final architecture:** keep one explicitly conditional coupled-action PDE as a mathematical test bed, preserve the distinction between real/complex sectors and kernel regimes, and stop before strict-source/physical promotion. The next dependent calculation is mechanical, not a new discovery: bind the real `A,B` and full complex Hessian, then evaluate the analytical semidefinite and Rayleigh tests and test Cauchy-domain admissibility. If no new typed nonlinear source map is provided, source verdict remains unchanged.
