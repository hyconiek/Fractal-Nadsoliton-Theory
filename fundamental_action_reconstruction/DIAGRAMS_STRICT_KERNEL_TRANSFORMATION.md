# Diagrams: legacy kernel -> strict kernel transformation

Status: `DRAFT_V0_RELEASE_SCAFFOLD__FINITE_BRIDGE_LEDGER_NOT_THEOREM_CLOSURE`

This file is the strict-kernel counterpart of `DIAGRAMS_KERNEL_TRANSFORMATION.md`. It is deliberately written as a release scaffold: it records the current finite construction of the bridge from the restored legacy intermediate kernel to the strict working kernel, while preserving the guardrails that block silent identity, silent role transfer, QW-2191 discharge, and ToE closure.

## 1. Ontological placement

The nadsoliton is treated here as the primordial information of the universe in a solitonic state. There is no deeper informational substrate under the nadsoliton. The preferred internal ordering remains:

```text
nadsoliton -> light -> matter -> emergent observer
```

Within that ontology, the kernel is not a separate entity below the nadsoliton. It is a compressed operational trace of nadsoliton structure. The legacy kernel was the first effective trace. The strict kernel is the later enriched trace after adding strict-side characteristics of the nadsoliton: nonlinear compression/damping, phase/topological data, and finite bookkeeping of carrier/node/edge information.

## 2. Historical construction line

```text
[legacy physical ansatz]
        |
        v
K_legacy_ont(d) = alpha_geo*cos(omega*d + phi)/(1 + beta_tors*d)
        |
        |  restored as intermediate incomplete bridge kernel
        |  not a final co-equal strict kernel
        v
[explicit completion ledger]
        |-- amplitude scalar normalization A(d)
        |-- phase/frequency transport P(d)
        |-- damping/compression passage D(d)
        |-- finite diagonal quotient Q = diag(K_strict/K_legacy)
        |-- APD symbolic cancellation certificate
        |-- GF(2) phase/topology certificates
        v
K_strict_gate(d) = cos(omega*d + phi)/(1 + beta*d^eta)
        |
        |  completed/enriched strict working kernel
        |  still missing source/selector/role-transfer theorems
        v
[frontier: strict sources + selector + role transfer]
```

The historical point is asymmetric. `K_legacy_ont` is restored as an intermediate bridge kernel because earlier work indicates that the strict kernel is an enriched completion of the legacy kernel. But the legacy kernel is not allowed to silently substitute for the strict kernel: it did not encode the full strict nadsoliton characteristics, especially nonlinear compression and certified GF(2)/cohomological phase data.

## 3. Kernel split and bridge equation

The audited split is:

```text
K_legacy_ont(d) = alpha_geo*cos(omega_L*d + phi_L)/(1 + beta_tors*d)
K_strict_gate(d) = cos(omega_S*d + phi_S)/(1 + beta*d^eta)
```

The finite bridge ledger packages the comparison on the audited `Z12` domain by factors

```text
A(d) = 1/alpha_geo
P(d) = cos(omega_S*d + phi_S)/cos(omega_L*d + phi_L)
D(d) = (1 + beta_tors*d)/(1 + beta*d^eta)
```

under explicit admissibility conditions `alpha_geo != 0`, `cos(omega_L*d + phi_L) != 0`, `1 + beta_tors*d != 0`, and `1 + beta*d^eta != 0` on the finite domain. The formal cancellation target is:

```text
K_legacy_ont(d) * A(d) * P(d) * D(d) = K_strict_gate(d)
```

This is a finite bridge-comparison certificate, not a derivation of the factors from strict nadsoliton dynamics.

## 4. What the strict completion adds over legacy

The strict kernel adds at least the following strict-side data, including nonlinear damping/compression as a strict-side addition:

| Strict addition | Current finite status | Why it is not silent legacy inheritance |
|---|---:|---|
| Amplitude scalar normalization | exact scalar witness `alpha_geo=4 ln(2)` as normalization, not role theorem | removes the legacy global scalar but does not transfer electroweak roles |
| Phase/frequency transport | affine transport certificate | aligns phase coordinates but does not derive an omega/phi source |
| Nonlinear damping/compression | `beta=1`, `eta=9/5` finite identifiability and linear-vs-nonlinear separation | legacy denominator `1+beta_tors*d` cannot equal `1+d^(9/5)` at multiple positive nodes with one constant gamma |
| Finite diagonal quotient | unique `Q=diag(K_strict/K_legacy)` on audited `Z12` support | a comparison map, not raw identity |
| GF(2) phase pattern | unique four-edge solution `1->2`, `5->6`, `7->8`, `9->10` | legacy linear denominator is sign-positive and cannot supply phase flips |
| Path/cycle cohomology | path `H^1=0`; closed cycle has one odd-parity class | the `d=0` anchor is a `C0` gauge fix, not the `C1/im(delta)` generator itself |

## 5. Computed bridge-status diagram

```text
finite bridge components
  amplitude scalar normalization          [certified finite witness]
  phase/frequency affine transport        [certified finite witness]
  damping/compression identifiability     [certified finite witness]
        |
        v
  APD assembly + symbolic cancellation    [certified finite comparison]
        |
        v
  strict kernel comparison on Z12          [available]
        |
        +--> strict source theorem         [open]
        +--> chi_11 selector/source        [open]
        +--> role-transfer theorem         [open]
        +--> ToE closure                   [open]
```

Thus the current honest claim is: the repository has a finite, auditable bridge-comparison ledger from the legacy intermediate kernel to the strict kernel. It does not yet have a strict dynamical source theorem for the bridge factors or the selector bit.

## 6. Role-transfer consequence

Because the strict kernel is richer than the legacy kernel, role transfer must be updated rather than copied. Once the legacy -> strict bridge is fully specified, the next audit must classify each legacy physical-role claim as one of:

1. **survives unchanged** under the strict completion;
2. **survives only in modified strict form** because nonlinear compression/topology changes the semantics;
3. **is rejected** as a legacy-only interpretation.

Until that audit is theorem-grade, the following legacy roles remain blocked from silent transfer: `sin^2(theta_W)=alpha_geo/12`, `alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)`, the `beta^N` gravity hierarchy, and the candidate `beta_tors -> chi_11` selector story.

## 7. Hard limits preserved by this scaffold

- No unqualified identity `K_legacy_ont == K_strict_gate` is claimed.
- No derivation of `A(d)`, `P(d)`, `D(d)`, `omega/phi`, `beta/eta`, or the transport cocycle from strict nadsoliton dynamics is claimed.
- No `beta_tors -> chi_11` theorem is claimed.
- No legacy physical-role transfer to `K_strict_gate` is claimed.
- No `QW-2191` selector discharge is claimed.
- No ToE closure is claimed.

## 8. Next release-build tasks

1. Promote this scaffold into a full release diagram file with references to every certificate report.
2. Add theorem-grade derivations or explicit axioms for the strict source atoms: APD source, phase/frequency source, damping beta/eta source, and `chi_11` selector source.
3. Complete the separate strict role-transfer audit after the bridge is no longer only finite-comparison scope.
4. Connect the strict Lagrangian/EOM draft to 4D covariant residual-zero witness tables.
