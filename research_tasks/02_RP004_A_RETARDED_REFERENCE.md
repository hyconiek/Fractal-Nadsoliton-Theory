# RP004-A — One finite-window retarded outward reference bound

Execution prompt for GPT-6.1 Sol. **UNEXECUTED; P1, optional field track**.
Cwd: `/media/magazyn/home/teresa/Dokumenty/HYC/TOE/edison`.
Write only `fin_rp004_flux_continuation_20261019/RP004-A/`.

## 1. Problem

For the already fixed RP-004 neutral datum, bound the signed flux of the
actual zero-data retarded scalar reference u on the ray R=T/6. Prove a
strict finite-window lower contribution or identify the exact obstruction.

## 2. Importance

Positive on-shell Gamma does not prove transported energy from a prepared
Cauchy problem. This is the missing linear-reference component of the
existing nonlinear flux gate, not a new radiation campaign or natural source.

## 3. Exact inputs

- `AGENTS.md`, `research_tasks/00_EXECUTION_RULES.md`.
- `research_tasks/RP004_SHARED_CONTRACT.json`.
- `FIN_RP004_RETARDED_FLUX_CONTINUATION_20261019.md`.
- `fin_coherent_mode_radiation_campaign_20261017/MODE_AND_FLUX_CONTRACT.json`.
- `fin_coherent_mode_radiation_campaign_20261017/ON_SHELL_CERTIFICATE_OR_FAILURE.md`.
- `fin_coherent_mode_radiation_campaign_20261017/OUTWARD_FLUX_THEOREM_OR_OPEN_LEMMA.md`.
- `fin_coherent_mode_radiation_campaign_20261017/PHASE_AND_INITIAL_SLIP_CONTROLS.md`.
- `fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md`.
- `fin_coherent_mode_radiation_review_20261019/RETARDED_FLUX_PROOF_ROUTE.md`.

## 4. Prior work and duplicate check

RP-002 bound mode and RP-003 scattering normalization, both roots and
Gamma>125/73728 are already paid in the restricted family. The 26-file
radiation search found the current target explicitly open; the old QW-467
numeric river-model flux is a different object, not its solution.
Before work, use targeted `rg -n 'OUTWARD_FLUX|retarded|limiting.amplitude'`
in the RP source, review and any existing continuation directory. Reuse a
later equivalent proof if present; do not repeat the Gamma calculation.

## 5. Paid equations

L=-d_x²+4-6sech²x, a in[-3/4,-1/4], Omega=2sqrt(a+2),
S=(9/4)tanh(x)sech^4(x), and

    u_tt+Lu=-S(1+cos Omega t), u(0)=u_t(0)=0.

u is odd. The odd shape mode is retained; the even translation mode is
absent by the actual parity. Continuum momentum q in[1,sqrt3] is away from
threshold. The Jost/Darboux resolution and Gamma normalization are inputs.

## 6. Missing statement

Target a bound on `F_u(R,T)=int_0^T[-u_tu_x](t,R)-[-u_tu_x](t,-R)dt`
with R=T/6, uniform in the frozen a interval. A sufficient outcome is
F_u>=c_ref T-E_ref(T), c_ref>0, E_ref=o(T), with the threshold and dependence
paid. This is a target, not an assumption or a promised theorem.

## 7. Method

Use the exact spectral Duhamel formula, retaining the DC term and the
finite-time removable singularity at energy Omega². Decompose into continuous
and discrete components only with their signed cross terms accounted for.
Try uniform oscillatory-integral/limiting-amplitude estimates inside the
group-velocity cone, with turn-on and exponentially localized source tails.
Group velocity q/Omega is not a universal sharp propagation front.

## 8. Actions

1. Write the retarded spectral solution and verify its sign/zero Cauchy data.
2. Isolate the continuum pole contribution with both directions and correct
   dq/(2pi) normalization.
3. Bound the finite-time transient on the moving boundary ray, including
   thresholds, DC, odd bound mode and source tails.
4. Produce one lower bound with all errors, or the first unpaid integral.
5. Export the exact norm/time requirements needed from RP004-B/C.

## 9. Assumptions

No eternal outgoing wave may replace u. No independent initial radiation,
parameter retuning, added Phi6 or modified kernel is allowed. Work under the
shared contract; active e1 is spectator-kernel-blind. This is classical
conditional mathematics on a supplied 1+1 domain.

## 10. Result types

Finite-time lower-bound theorem; counterexample to the frozen reference
target; or exact open oscillatory-integral bound. A fixed-R long-time result
alone is partial if its constants do not pay the ray R=T/6.

## 11. Acceptance and independent checks

Verify the reference local energy identity independently of the Jost formula.
Check signs/factors against the inherited explicit spectral resolution.
Numerical quadrature may falsify a mistaken formula, never certify uniform
ray asymptotics. Keep the circular-mode/DC distinction and initial-slip control.

## 12. Dependencies

No result from ID-001 or RP004-B is needed for this reference problem.
A and B may run in parallel after selecting the field track, with disjoint
writes. Its output is not an exact nonlinear flux theorem.

## 13. Failure

Return the failed ray/transient bound. Do not choose a numerically favorable
window. One analytically justified alternative window may be proposed for
architect review, not silently installed into B/C.

## 14. Stop

One spectral route and one concrete repair. No Gamma precision work, broad
parameter scans, full scattering programme or absorption by a second wall.

## 15. Deliverables

`RETARDED_REFERENCE_BOUND_OR_OPEN_LEMMA.md`, `REFERENCE_ERROR.json`,
`NORMALIZATION_AND_RAY_CONTRACT.json`, common proof/report/handoff files.
Verdict: `RETARDED_REFERENCE_PAID`, `PARTIAL_FIXED_WINDOW_ONLY`,
`COUNTEREXAMPLE` or `OPEN_LEMMA`. Return; do not activate C without B.
