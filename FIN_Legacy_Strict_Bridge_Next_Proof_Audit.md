# FIN Legacy–Strict Bridge Next Proof Audit

## Metadata

- Date UTC: `2026-07-19T19:00:58.050296+00:00`
- Git HEAD at generation: `b7d112f091aff320d26b722810dcfc1ccdbea8fd`
- Python: `3.14.4`
- NumPy: `2.5.1`

## Executive verdict

This is a narrower proof-computational follow-up to the legacy-star audit. It constructs the explicit scalar/envelope/phase quotient objects between the reconstructed `K*` legacy class and `K_strict_gate`, then separates finite algebraic witnesses from missing source theorems.

- **[Proven finite obstruction]** Positive scalar and positive damping/envelope maps cannot fix the `K*` to strict sign mismatch; a phase/frequency transport object is mathematically necessary.
- **[Proven quotient fact]** `A` is a scalar multiplier: it drops out of normalized shape/sign comparisons, but it is still needed for unnormalized coupling strength, legacy `alpha_geo` roles, moment magnitudes, and any future action-density normalization.
- **[Strong evidence]** `K*` is easier to bridge than the rejected double-damped product because it is not exponentially dead and admits a clean positive damping quotient. It is still not close to strict without a sign-changing phase quotient.
- **[No-go on current artifacts]** The constructed quotient `Q(d)` is an algebraic finite assembly witness, not a source theorem, selector theorem, role-transfer theorem, unit theorem, or ToE closure.

## Prior repo anchors used

| Anchor | Imported status |
| --- | --- |
| P2760 | seven open gaps: ontology-to-kernel measure, amplitude normalization, phase/frequency/topological source, damping/compression bridge, moments-to-couplings provenance, Lagrangian reverse closure, stale flags |
| P2691 | alpha_geo=4 ln 2 and scalar-shape normalization real; no role-safe amplitude absorption source |
| P2692 | positive beta orbit/gauge representability real; no target-independent positive beta/Z_beta source theorem |
| P3167_P3168 | no strict unit source S+ on current monomial/state map artifacts |
| scratch_bridge_assembly | existing finite assembly quotient maps legacy to strict on audited nodes but explicitly does not export source, selector, role-transfer, or ToE closure |

## Kernel definitions compared

```text
K*_legacy(d;A,beta_l) = A*cos(pi*d/4 + pi/6)/(1+beta_l*d)
K_strict_gate(d)      = cos(0.18575*d + 0.16250)/(1+d^1.8)
K_rej(d)              = exp(-2.9*d)*(1+0.2*d)/(1+d)*cos(pi*d/4+pi/6)
```

## Sign obstruction theorem

For every A>0 and every positive radial envelope E(d)>0, sign(A*E(d)*cos(pi*d/4+pi/6)) = sign(cos(pi*d/4+pi/6)); therefore positive amplitude/envelope completion cannot map the K* sign row to the strict sign row.

| Row | Signs d=1..12 |
| --- | --- |
| strict | [1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1] |
| K* | [1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1] |
| K_rej | [1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1] |
| K* mismatch distances | [2, 3, 4, 5, 8, 9] |

**Consequence:** the bridge cannot be only `A` plus a positive damping/compression map. The missing theoretical object is a phase/frequency transport/source, not another scalar amplitude fit.

## Damping/compression quotient object

For `beta_l=0.05`, define `D_beta(d)=(1+beta_legacy*d)/(1+d^1.8)`. On `d=1..12` it is positive: `True`, monotone decreasing: `True`, with `D(1)=0.525` and `D(12)=0.0180578`.

This is useful because it cleanly separates linear legacy damping from strict nonlinear compression after phase and amplitude have been factored away. It does **not** derive `eta=1.8`; it only says what positive envelope quotient would be needed if the phase/source layer were supplied.

## Finite assembly quotient

After separating amplitude, phase, and damping, the exact finite quotient is:

```text
Q(d) = A^-1 * [cos(theta_strict(d))/cos(theta_legacy(d))] * [(1+beta_l*d)/(1+d^1.8)]
K*_legacy(d;A,beta_l) * Q(d) = K_strict_gate(d)
```

For `A=2.9`, `beta_l=0.05`, the maximum absolute identity error on `d=1..12` is `5.551e-17`. This is an exact finite algebraic witness, but the phase factor ranges from `-1.78262` to `3.63177` and has `4` sign changes, so it is not a positive dissipative bridge.

## What is `A`?

A is a scalar coupling/measure-density normalization multiplying the whole legacy profile; it is not a shape parameter.

| Question | Answer |
| --- | --- |
| Does A drop from normalized shape? | True |
| If A=1, is it needed? | If one works only in quotient-normalized shape space, A=1 is a valid gauge. If one claims alpha_geo roles or absolute physical couplings, A cannot be deleted; it becomes an open source/normalization obligation. |
| A is needed for | ['absolute unnormalized coupling strength', 'legacy role formulas that use alpha_geo explicitly', 'kernel moments before quotient/row normalization', 'any future action-density normalization theorem'] |
| A is not needed for | ['sign pattern', 'normalized profile shape', 'correlation/cosine similarity', 'row-stochastic Markov operator after full row normalization'] |

The value `A=4 ln 2=2.772588722239781` remains **[Proven]** as `ln 16`, the entropy of a uniform four-bit/16-state cell, but **[No-go]** as a physical unit, selector, or role-transfer source on current artifacts.

## Is the new legacy kernel easier to bridge to strict?

| Comparison | relative L2 | cosine | scalar |
| --- | --- | --- | --- |
| unconstrained scalar K_rej -> strict | 0.491662 | 0.870786 | 52.6998 |
| positive scalar K_rej -> strict | 0.491662 | 0.870786 | 52.6998 |
| unconstrained scalar K*(beta=0.05) -> strict | 0.997242 | -0.0742191 | -0.0208758 |
| positive scalar K*(beta=0.05) -> strict | 1 | -0.0742191 | 0 |
| best grid beta K* -> strict | 0.997201 | -0.0747643 | beta=0.0741925, A_fit=-0.0233015 |

The scalar table is intentionally harsh: because the strict row is positive until d=7 while K* has negative entries at d=2..5, a positive scalar-only fit to K* collapses to zero. The rejected product can be inflated by a huge scalar to match mostly the first entry, but this is not a bridge; its norm is tiny and it is exponentially dead for d>=2. K* is therefore easier only after allowing the necessary phase/frequency object O_P plus the clean damping quotient O_D.

K* is easier than K_rej for bridge work because it preserves a non-dead long-range hyperbolic profile and separates damping from phase; it is not close to strict without a nontrivial phase/frequency map.

## Constructed theoretical-object matrix

| Object | Constructed here | Status | Missing for closure |
| --- | --- | --- | --- |
| O_A scalar amplitude/measure normalization | quotient factor A^{-1}; proof that A drops from normalized shape | witness/gauge decomposition only | source theorem fixing A as alpha_geo or proving amplitude absorption is role-safe |
| O_D damping/compression map | D_beta(d)=(1+beta_l*d)/(1+d^1.8), positive and decreasing on d=1..12 | finite/exact envelope witness after phase separation | dynamical source of eta=1.8 and target-independent beta/Z_beta theorem |
| O_P phase/frequency transport | P(d)=cos(theta_s)/cos(theta_l) finite quotient on d=1..12 | algebraic quotient with sign changes, not source law | phase/frequency/topological source localizer and selector-safe transport theorem |
| O_S selector/orientation source | positive scalar/envelope sign no-go matrix | obstruction: QW-2191 remains open | non-premise orientation/polarity source plus coupling-polarity theorem |
| O_U unit/action-density bridge | none; imported no-go status from guardrails/P2760-P2764/P3167-P3168 | no-go on current artifacts | canonical reference-cell/action-density theorem with physical units |

## Recommended next honest step

**Recommendation:** Do not run another generic bridge replay. Attack exactly O_P: construct or refute a phase/frequency transport source localizer for P(d)=cos(theta_s)/cos(theta_l), while conditioning O_A and O_D as quotient witnesses only.

**Why:** The new K* already repairs the historical envelope, so the remaining hard obstruction is not scalar A and not positive damping alone; it is the sign-changing phase/frequency/source layer plus the still-open unit/selector obligations.

**Success criterion:** A machine-checkable theorem or finite exhaustive certificate that P(d) comes from a strict exported source object and is selector-safe; otherwise emit a bounded no-go for this exact P object.

Do not promote this audit to strict bridge closure, selector closure, role transfer, physical units, `L_total`, or ToE.
