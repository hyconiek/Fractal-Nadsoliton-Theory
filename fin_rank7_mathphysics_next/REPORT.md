# FIN MP7 mathematical-physics campaign — final scientific report

Date: 2026-09-21  
Namespace: **MP7-001–MP7-048**  
Working predecessor: `FIN_R7O3_TARGETP_HANDOFF_20260920.zip`  
Predecessor SHA-256: `7ed1f1de9640c1bbc353b1360cafb90180c6772bf6ec8920f358906e07ab22ad`

## 1. Strongest new conclusions

The strongest result of this campaign is not another coverage percentage. It is an analytic and interval-assisted reduction of the full minimization problem followed by a global transition theorem.

First, MP7-007–011 prove an all-orders positive Fourier expansion of the exact partition function and the phase-alignment inequality `Z(phi;b)<=Z(0;b)`, with equality classified on every zero-amplitude stratum. Consequently **every global minimizer of the full seven-dimensional dual is D12-equivalent to an aligned nonnegative four-amplitude representative**. This is a nonlinear partition-function theorem; it does not assert that every stationary point is aligned.

Second, MP7-016 exhausts the aligned stationary problem at the exact gain `g=37/10`. There are exactly three aligned nonnegative roots: the uniform root, one localized minimum and one index-one saddle. The nonzero energies satisfy

- localized: `[0.008235102001579031, 0.008235181373454593]`,
- saddle: `[0.0487900472912212, 0.04879010290166141]`,

while the uniform energy is exactly zero. Combined with phase alignment, **the uniform state is the unique full-X7 global minimizer at g=3.7**.

Third, MP7-017 globally exhausts the certified equal-energy event. The event lies in

`g_eq in ['3.7183448971203875136738542882050261929327867539951469340619588641', '3.7183448991203875136738542882050261929327867539951469340619590691']`.

At the event only the uniform state and the 12 translated localized minima attain the global value; the saddle stays above `0.046555431928778354`. The primal identity `V_g(p)=V_geq(p)+(g_eq-g)||X7^T p||^2/2` then proves that the uniform state is uniquely global for every smaller gain. Thus, **for each supplied strict spectral tuple, its certified localized/uniform equal-energy event is its first global transition in this finite model**.

## 2. Target-P curvature result and normalized margin

The R7O3 Target-P candidate was intake-audited rather than regenerated. MP7 uses decreasing covariance eigenvalue order `lambda1>=lambda2>=lambda3>=lambda4`; the invariant statement is that at most one eigenvalue exceeds `tau0=67/250`.

MP7-020 converts the heterogeneous compact certificates and accepted tails into a common normalized spectral statement. Over all **25656** compact terminals plus the accepted unbounded tails,

`lambda2(M4) <= 0.267999946371058`,

with the rigorous positive separation

`m >= 6500000000000000000000000000/121203210760018863485060548125326941 ~= 5.36289423295e-08`.

The weakest compact cell is an R7O3 repair. The accepted tails have much larger separation (`~0.000556756`), so they do not determine the global minimum margin. This also resolves a previous interpretation issue: the old smallest PD-test quantity was a witness/minor margin, **not** itself a covariance spectral gap.

MP7-019 proves the correct invariant formulation under linear coordinate changes: covariance and Hessian matrices transform by congruence together with the quadratic metric, so generalized eigenvalues/inertia are the meaningful quantities. The threshold `67/250` is tied to the supplied feature normalization and must not be advertised as a normalization-free physical constant.

## 3. Boundary supports and the stationary 4D-to-7D bridge

MP7-013 identifies the only structurally possible aligned supports; MP7-014 excludes every nonzero incomplete support for `0<g<=250/67`. Hence any nonzero aligned stationary point in this window is interior.

MP7-037 proves analytically that all three odd/angular Hessian directions are positive at an interior aligned stationary point. MP7-038 combines this with Target P to obtain a scoped full-X7 result: **the full mediator Hessian has negative index at most one at any interior aligned stationary point for `0<g<=250/67`**.

This is deliberately not an unrestricted stationary theorem. MP7-039 studies the known two-harmonic counterexample family instead of erasing it. A unique local transverse crossing is certified at

`g* in [5.1718418319285595, 5.171841831957076]`

with derivative `[0.0016957351021672582, 0.0016960387698923191]`. Two symmetry-related blocks cross simultaneously, giving a two-dimensional critical space carrying the standard `D3` representation; `Re(z^3)` is symmetry-allowed. Therefore the event is not a generic one-dimensional pitchfork.

## 4. Fold, response and conditional kinetics

MP7-025/026 turn the earlier simple-fold certificate into controlled local asymptotics. With `epsilon=g-g_fold`, for

`0 < epsilon <= 1e-7`,

MP7-026 proves exactly two local reduced roots and obtains

- `|xi|/sqrt(epsilon)` in `[1.8619243247549937, 1.9186338692641498]`,
- local saddle-minus-minimum `Delta Phi / epsilon^(3/2)` in `[0.5126867995891401, 0.558650696110349]`,
- `|lambda_soft|/sqrt(epsilon)` in `[0.21490520515383432, 0.2349154318867044]`.

The local energy difference is not claimed to be a global escape barrier.

MP7-022/023 separate source conventions and quantify the collective mode. At coexistence the four components of `ds/dg` are all positive with the narrow intervals stored in `results/MP7-022_023_quantitative_response.json`. At the fold the leading covariance eigenvalue is separated from the global upper bound on the second eigenvalue by at least

`0.016442902010540756517838171451187902090481506254571150414122865450852446446645661`.

MP7-028/029 then add explicit gradient dynamics. Hessian inertia fixes local instability count but mobility fixes rates. For `L=I` the stable-fold relaxation time grows like `epsilon^(-1/2)` with controlled prefactor; changing `L` changes the rate without changing the equilibria. This is a demonstration of kinetic nonuniqueness, not a derivation of a physical clock.

## 5. Finite-N statistical extension

MP7-031–033 define an explicit, added finite-copy model. They prove the exact occupation-number Gibbs weight, the variational limit to `V_g`, the seven-dimensional Gaussian auxiliary representation, the probability/mediator determinant identity and exact finite-N fluctuation-response relations. MP7-036 checks these identities independently at small `N` without Monte Carlo; every validation gate passes.

At coexistence, MP7-034 proves that the localized minimum has a D12 orbit of size 12. The total Gaussian prefactor ratio of the localized family to the uniform minimum at equal energy is

`[0.5437337926269059, 0.5437338754161343]`.

Because this ratio is below one, equal **finite-N mass** occurs at a slightly larger gain than equal energy. MP7-035 gives

`g_N = g_eq + c*/N + o(1/N)`,

with

`c* in [1.3467658209548814, 1.3467661619763278]`.

After MP7-017 this is a global leading equilibrium asymptotic, not merely a comparison of two arbitrary local roots. The campaign also supplies explicit fixed-`g_eq` local Laplace-error thresholds. For both caps together, sufficient `N` values are approximately:

| Relative local error | sufficient N |
|---:|---:|
| 25% | 6.08412e+08 |
| 10% | 4.41389e+09 |
| 5% | 2.05878e+10 |
| 1% | 7.55279e+11 |

These huge numbers are conservative sufficient bounds produced by crude uniform inequalities. They are **not** evidence that the numerical asymptotic only becomes useful at those N. The explicit global outside-cap remainder uniform in `g=g_eq+c/N` remains the principal open quantitative finite-N atom.

## 6. Phase-census scope

The accepted frozen-amplitude fixture has exactly 60 phase critical points. MP7-040 checks what those roots actually mean: none satisfies full positive-g stationarity at the same frozen amplitudes. Thus they are a phase census, not 60 full equilibria or physical states.

MP7-041 then chooses one symmetry-locked root and proves an exact constant phase continuation on a ±0.1% amplitude box. This converts the nearby full-equilibrium search on that family into a four-amplitude radial problem; it does not prove a nearby equilibrium exists or is globally minimizing.

## 7. Robustness

MP7-021 distinguishes trivial scaling equivalence from genuine coefficient perturbation. With unscaled fields held fixed, the global Target-P margin gives the conservative sufficient relative spectral-weight increase

`<= 2.0010803381e-07`.

At the localized coexistence minimum, the local differential sensitivity gives a finite operator bound recorded in `results/MP7-021_local_spectral_sensitivity.json`. This is a local derivative bound and must not be read as a large finite perturbation radius.

## 8. Reproducibility and negative controls

MP7-043 records an acyclic proposition dependency graph. MP7-044 exercises the required semantic and implementation failures, including metric omission, witness rescaling, phase-root/full-root confusion and M4/M7 substitution.

MP7-045 runs the new executable results from a fresh copied input tree with parameterized roots. **19 JSON outputs were compared and all scientific fields matched.** The only ignored difference is the wall-clock `elapsed_seconds` field in MP7-020. This establishes clean-environment reproducibility of the same implementations; it is not claimed as a fully independent second implementation.

## 9. Scientific status

The campaign therefore produces several new analytic and interval-assisted results, including phase alignment, exact global minimization at `g=3.7`, the first global transition, a normalized global Target-P margin, controlled fold scaling, a scoped stationary 4D-to-7D bridge, and an explicit finite-N equilibrium extension.

The result register is `CLAIM_REGISTER.json`. Explicit nonconclusions are in `NONCONCLUSIONS.md`. The finite follow-up list is `NEXT_ATOMS.md`.

The main remaining mathematical atom is **not** a hidden condition on the global-transition theorem. It is the stronger finite-N quantitative problem of a global outside-cap error bound uniform in the `g_eq+c/N` scaling window.
