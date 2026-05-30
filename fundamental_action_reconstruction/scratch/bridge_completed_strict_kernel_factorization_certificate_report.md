# Completed strict-kernel factorization certificate

Status: `strict-kernel-current-full-form-with-explicit-legacy-completion-factorization`

## Kernel reading

- Current live kernel: K_strict_gate(d)=cos(omega_S*d+phi_S)/(1+beta*d^eta)
- Legacy reading: K_legacy_ont is the historical/nadsoliton-characteristic carrier, not the live current kernel.
- Completion reading: K_strict_gate is represented as legacy carrier plus explicit completion factors for alpha removal, phase/frequency transport, damping/compression, and the still-open orientation bit lane.
- Guarded identity: Pointwise factorization is exact where cos(omega_L*d+phi_L) != 0; it is not yet an internal theorem deriving the factors from nadsoliton dynamics.

## Factorization summary

- Domain d values: `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]`
- Max abs residual from legacy full: `1.110e-16`
- Max abs residual from legacy reduced: `1.388e-17`
- Residual tolerance pass: `True`
- Damping compression positive/decreasing: `True` / `True`
- Damping factor d=1 to d=11: `0.505000000000` -> `0.014623674672`
- Damping factor ratio d1/d11: `34.533044`
- Phase transport sign pattern: `[1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`

## Completion factor interpretation

- `amplitude`: 1/alpha_geo removes the legacy explicit amplitude; role-safe use still requires bridge discipline.
- `phase_frequency`: cos(strict phase)/cos(legacy phase) is the exact pointwise transport from legacy torsion/resonance phase to strict gate phase on the sampled domain.
- `damping_compression`: ((1+beta_tors*d)/(1+beta*d^eta)) is the exact pointwise compression factor upgrading legacy hyperbolic damping into strict eta damping.
- `orientation_bit`: The factorization does not supply beta_tors -> chi_11; that remains the one-bit theorem-target from the reactivation audit.

## Proof certificate

- `factor_identity`: For every listed d, K_strict(d) - C_full(d)*K_legacy_full(d) = 0 up to floating arithmetic, and similarly for the reduced no-alpha carrier.
- `completion_objects`: C_full=(1/alpha_geo)*phase_transport*damping_compression and C_reduced=phase_transport*damping_compression are explicit completion factors.
- `computational_positive`: The current full kernel can be represented exactly as strict completion of the legacy carrier on the finite d=1..11 audit domain.
- `theoretical_limit`: This is an exact factor certificate, not a derivation of the factors, not the beta_tors->chi_11 theorem, and not QW-2191 closure.

## Hard limits

- No claim that the legacy kernel is the current live kernel; K_strict_gate is the current full form.
- No unqualified identity K_legacy_ont == K_strict_gate is asserted; only explicit factorization/completion is certified.
- No beta_tors -> chi_11 theorem is asserted.
- No proof derives the completion factors from strict nadsoliton dynamics yet.
- No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
