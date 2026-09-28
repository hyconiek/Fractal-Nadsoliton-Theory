# FIN PHYSICAL BRIDGE CAMPAIGN I — HANDOFF AFTER PHYS-001–007

Date: 2026-09-29
Baseline repository commit: `97c4231f33632800fd817fe2294555bb8bcb041f`
Scope: PHYS-001–007 only. **PHYS-008 was not executed.**

## Executive verdict

Campaign I closes a rigorous *model-class bridge*, not a fundamental-physics bridge.

1. **PHYS-001 — EXACT.** The finite-N Gibbs core is exactly a generalized 12-state matrix Curie–Weiss model. Because `A7=X7 X7^T`, it is equivalently a discrete vector-spin mean-field model with 12 allowed internal vectors in `R^7`. This is not seven-dimensional physical space.
2. **PHYS-002 — scoped PASS with source gap.** A fail-closed small-N implementation reconstructed from versioned sources passes A7, generator, detailed-balance, C12-sector and fixture tests. `rho_3=0.13143978619564556`, differing from the declared fixture by `2.22e-16`. The historical untracked `FIN son` files could not be byte-for-byte audited; no claim of repairing those originals is made.
3. **PHYS-003 — PASS.** The exact static fingerprint is `S_k(0)=1`, `dS_k/dg|0=(N-1) Lambda_k/(12N)`. FIN has zero first-order sectors `k=1,2` and active `k=3,4,5,6`; after scale normalization the active shape ratios are `1.12142406146`, `1.17191711554`, `1.19413370400` relative to `Lambda_3`. N=2 additionally gives an exact finite-g pair histogram, so the first physical test need not rely on a weak-g approximation.
4. **PHYS-004 — PASS with a negative dynamic conclusion.** Trace-matched full Potts, flat P7, three active-weight perturbations and a k1/k2-leakage control were frozen. Static pair distributions distinguish the matrices. Dynamic rates do **not** belong to A7 alone: at the same A7 and g, heat-bath, Metropolis and Barker preserve the same Gibbs law but have materially different sector rates after a common g=0 clock normalization.
5. **PHYS-005 — NO_NEW_SOURCE.** Current inputs supply the strict-kernel parameters and the active support selection. D12/PSD/rank/trace do not uniquely determine the three remaining active spectral ratios. A Gaussian-mediator parent realizes any suitable PSD matrix, but is only an engineered realization unless `K` and `c_j` are independently fixed/measured before seeing the FIN fingerprint.
6. **PHYS-006 — DESIGN_ONLY scoped feasibility.** A programmable electronic/mixed-signal categorical sampler with explicit controller and physical stochastic source is the cleanest realization contract. Existing multi-state p-bit work supports feasibility of stochastic multi-state/Boltzmann hardware, not an already demonstrated exact q=12 A7 leave-one-out FIN device. Ring-oscillator Potts optimization is not enough to establish Gibbs sampling/rates.
7. **PHYS-007 — conditional preregistration PASS.** Freeze the first test as an N=2 conditional pair-difference histogram at `g=3`, with `g=0` negative control. Primary alternatives are full Potts, flat P7 and two 10% trace-preserving active-spectrum perturbations. The minimum ideal primary TV separation is `0.01178719`. Requiring a pre-validation total calibration envelope `TV<=0.003` per model leaves positive robust separation `>=0.00578719`. Closer 5%/2% alternatives remain sensitivity-only rather than being hidden by retuning.

## What is standard statistical physics

- Gibbs measure and multinomial degeneracy for finite mean-field copies.
- Generalized q-state Curie–Weiss / matrix-interaction mean-field models.
- Reversible heat-bath, Metropolis and Barker dynamics sharing a stationary Gibbs law while differing dynamically.
- Static fluctuation/response identities such as differentiation of Gibbs expectations by covariance.
- PSD factorization and Gaussian completion-of-square as a generic construction.

These facts do not select FIN.

## What genuinely depends on FIN/A7

- q=12 in the declared finite-N model.
- The Fourier cutoff/support `k=3,4,5,6` of the rank-seven construction.
- The concrete active eigenvalues and their three scale-free ratios.
- Consequently, the exact N=2 pair-difference distribution and the small-g `S_k` slope pattern.
- Any dynamic number obtained only after additionally naming a kinetic convention is FIN+A7+kinetics, not A7 alone.

## What has independent physical justification

At present, only the *class-level machinery* has independent physical justification: generalized spin/Curie–Weiss statistical mechanics is established physics, and published hardware demonstrates that multi-state stochastic/Potts-like samplers can be physically constructed. That supports feasibility of an engineered realization.

There is **no independent physical source currently demonstrated for the specific A7 kernel/support/eigenvalue tuple**. `G_FROZEN` has an internal operational origin (the documented barrier-balance construction), but that does not make it a constant of nature.

## What remains a manual/modeling assumption

- numerical strict-kernel parameters and the physical reason for them;
- selection of the active Fourier support `k=3,4,5,6`;
- the concrete three active spectral ratios after trace normalization;
- physical meaning of the 12 labels and of the seven feature coordinates;
- dimensionless gain g and its mapping to a real apparatus;
- update clock and choice of heat-bath versus other reversible kinetics;
- hardware realization of 12 arbitrary logits, physical RNG quality, delays, schedule and readout model;
- any claim that the programmed Hamiltonian equals actual thermodynamic/electrical energy of a device.

## Boundaries preserved

- Task 335 was not opened.
- Protocol 333 was not changed: `theta=2`, `Tprep=4`, deep seed, unrestricted preparation, no hard wall, no postselection.
- Task 337 remains WIP; no cross-N/cross-g response law was claimed.
- The current scoped repository intake still does **not** promote `Gamma=B4` to repository-certified status; the exact-input, LP and scalar-minimum enclosure obligations remain.
- No N>=4 production calculation was required; the campaign used N=1,2,3 only.
- No hardware purchase, lab contact, private data acquisition or AGENTS.md modification occurred.

## Architect checkpoint input

PHYS-008 should decide whether to proceed as an **engineered statistical-physics bridge** despite `NO_NEW_SOURCE`, or stop the fundamental-source lane until a genuinely independent A7 source atom exists. The execution agent stops here as instructed.
