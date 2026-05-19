# P2025 S975 Strict Cutkosky Same-Scheme CohomologyAmplitudeBridge Seed Theorem

Status: `OPEN_OBSTRUCTION_WITH_TRACE`.

v37 adds joint-fit perturbation robustness diagnostics (feature-matrix jitter),
exporting bounded perturbation residual span.

v37 also preserves the v36 joint-fit multistart robustness diagnostics, exporting
`multistart_residual_l2_span` as an additional bounded stability gate.

v36 also preserves the v35 channel holdout-rotation diagnostics for the joint coupled bridge fit,
exporting held-out channel residual bounds as an additional robustness gate.

v35 also preserves the v34 solver crosscheck
(`L-BFGS-B` vs `SLSQP`) and a `lambda_spread` sweep, exporting bounded
objective-gap and residual-span diagnostics.

v34 also preserves the v33 joint coupled bridge fit precursor for task-2/task-7
with an explicit spread-penalized objective solved by L-BFGS-B, exporting
joint residual and coefficient-spread diagnostics under OPEN status.

v33 also preserves the v32 cross-channel coefficient-spread gate for the task-2/task-7 bridge
fit, bounding inter-channel drift of shared-basis coefficients.

v32 also preserves the v31 unified strict-lane stability envelope for the task-2/task-7 bridge,
combining LOOCV, bootstrap-over-channels, and tolerance-sweep maxima into one
bounded gate (`phase_bridge_stability_envelope_max`).

v31 also preserves the v30 channelwise quadrature-tolerance sweep diagnostics for the strict-lane
phase-space Cutkosky precursor (task #2), exporting per-channel tolerance spans
and a bounded global tolerance span gate.

v30 also preserves the v29 bootstrap-over-channels stability diagnostics for the strict-lane
phase/common-basis link precursor, exporting bounded `bootstrap_residual_l2_p95`
under explicit non-closure status.

v29 also preserves the v28 leave-one-channel-out stability diagnostics for the strict-lane
phase/common-basis link precursor (task #2 -> #7 bridge step), exporting
held-out residuals and a bounded `loocv_residual_l2_max` gate.

v28 also preserves the v27 strict-lane phase/common-basis link precursor:
it projects channel phase-space aggregates onto a shared basis and exports
matrix coefficients, conditioning, and residual diagnostics under explicit
`OPEN_PRECURSOR_NOT_CLOSURE` status.

v27 also preserves the v26 strict-lane channelwise phase-space Cutkosky precursor for task #2:
it exports channel-parameterized phase-space integral tables over a shared
`s` grid and enforces positivity/monotonicity gates per channel.

v26 also preserves the v25 strict-lane DiscM common-basis integration precursor for task #7:
it exports a shared basis matrix, per-channel least-squares coefficients,
residual norms, and bootstrap uncertainty bounds, explicitly marked as
`OPEN_PRECURSOR_NOT_CLOSURE`.

v25 also preserves the v24 strict-lane QW-2191 selector-premise precursor for task #6:
it exports an explicit symmetry-breaking selector-weight ansatz scan and is
hard-marked as `NON_STRICT_PREMISE_PACKET`, with no strict-core closure claim.

v24 also preserves the v23 strict-lane PO2 sufficiency trace precursor for task #5:
it exports an explicit symbolic `L_total`, Euler-Lagrange traces,
`DELTA_BG_Yf` symbolic reduction under `C1..C4` constraints, Hessian-rank
certificate, and finite numeric replay rows confirming near-zero constrained
residual.

v23 also preserves the v22 strict-lane background-transport precursor for task #3
(`nu` branch): it exports symbolic FRW→Bianchi and Bianchi→FRW transport
maps, evaluates finite-grid bijective closure errors, and records a bounded
symmetry-commutator proxy over the `nu` scan.

v22 also preserves the v21 strict-lane constructive PO3 nonempty-certifier precursor (task #4):
it runs bounded L-BFGS-B over strict-kernel parameter space and exports
an explicit feasible candidate (`omega, phi, beta, eta`) together with
constraint checks (positivity floor, monotone phase profile, beta>0, eta>=1)
and a symbolic covariant-consistency proxy value.

v21 also preserves the v20 strict-lane backend renormalization-B1 precursor for task #1:
it builds strict-kernel B1 invariant features (`R2`, `Ric2`, `Riem2`, `GB`),
fits provisional backend coefficients (`a_R2`, `a_Ric2`, `a_Riem2`, `a_GB`)
with SciPy least-squares, and crosschecks them with a SymPy normal-equation
solve. The packet exports both coefficient sets, solver-gap, residual norms,
and the symbolic counterterm operator object.

v20 also preserves the v19 strict-lane multi-channel backend `DiscM_loop` precursor layer
(3 channels: graviton→gauge-gauge, graviton→fermion-fermion,
graviton→scalar-scalar) with per-channel dual-solver crosscheck
(Nelder-Mead + L-BFGS-B), explicit solver-gap bounds, and global nonzero
loss spread to prevent false closure.

All prior protections are preserved: 7-task ToE gap ledger,
explicit OPEN statuses, digest self-consistency, reproducibility locks,
FD-vs-symbolic checks, grid/tolerance robustness, bootstrap/conditioning,
holdout checks, and pre-closure nonzero-loss gating.

No `DiscM = CutSum` closure is claimed.
