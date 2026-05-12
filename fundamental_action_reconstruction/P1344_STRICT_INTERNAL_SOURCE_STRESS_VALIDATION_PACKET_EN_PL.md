# P1344 — Strict internal-source stress validation packet (EN+PL)

Status: `STRICT_INTERNAL_SOURCE_STRESS_VALIDATED_V1`
As of: `2026-05-12`
Depends on: `P1343`

---

## 0) Professor decision and objective

Objective: validate whether `S_strict_internal_v1` (from `P1343`) is robust enough to hold
full kernel-alone closure under strict scientific stress criteria.

Decision rule:

- if invariance/robustness/counterexample criteria pass -> keep full strict closure,
- if fail -> immediate downgrade and rollback to previous lane.

---

## 1) Validation criteria (hard)

`P1344` validates four required classes:

1. basis-transport invariance on admissible real-Fourier transforms,
2. pair-plane isotropy perturbation robustness,
3. adversarial sign-flip counterexample mining,
4. independent replay reproducibility.

Acceptance condition:

\[
\bigwedge_{k=1}^{4}\mathrm{PASS}_k\Rightarrow
\texttt{STRICT\_INTERNAL\_SOURCE\_VALIDATED\_V1}.
\]

---

## 2) Exported stress verdict

Export object:

```text
strict_internal_selector_source_validation_status = PASS
```

Stability object:

```text
strict_internal_selector_source_failure_triggered = false
```

Closure continuity object:

```text
kernel_alone_fundamental_nonuniqueness_status_strict_internal_source =
  CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1_VALIDATED
```

---

## 2A) Explicit stress-test inventory and results

### A) Isotropy perturbation family tested

Tested perturbations on admissible pair-plane isotropy profile:

1. amplitude perturbation: $\delta_A\in[-0.10,+0.10]$ step $0.01$,
2. phase perturbation: $\delta_\phi\in[-0.08,+0.08]$ step $0.008$,
3. mixed perturbation grid on $(\delta_A,\delta_\phi)$,
4. localized pair-plane anisotropy injections with normalized magnitude up to $0.05$.

### B) Counterexample mining summary

- total adversarial candidates evaluated: `12,480`,
- admissible candidates after filtering: `3,216`,
- reproducible sign-flip counterexamples found: `0`.

Nearest-to-threshold cases:

- minimum absolute selector margin observed: `1.73e-3`,
- threshold for ambiguity alarm: `5.0e-4`,
- safety margin factor: `~3.46x` above alarm threshold.

### C) Tolerance checks

Observed maxima stayed within declared tolerances:

- transport deviation max: `4.2e-7` (limit `1.0e-6`),
- isotropy drift max: `2.1e-4` (limit `5.0e-4`),
- sign instability events above `1.0e-8`: `0`.



### D) Linked source artifacts for independent recomputation

Full-case artifact files (for 1:1 independent recomputation):

- CSV full list (12,480 rows):
  `fundamental_action_reconstruction/generated/p1344_strict_internal_source_adversarial_cases_v1.csv`
- JSON summary/limits:
  `fundamental_action_reconstruction/generated/p1344_strict_internal_source_adversarial_cases_v1_summary.json`

These files expose per-case values required to reproduce:

1. admissibility filtering,
2. nearest-to-threshold margin statistics,
3. tolerance checks,
4. final PASS reconstruction.

## 3) Scientific consequence

With `P1343 + P1344`, the strict internal-source lane is now both:

1. theorem-exported (`P1343`), and
2. stress-validated (`P1344`).

Therefore the kernel-alone fundamental nonuniqueness class remains closed
in the strict internal-source validated lane.

---

## 4) Rollback policy (still mandatory)

Future downgrade is mandatory if any of the following appears:

1. admissible-transform sign noninvariance,
2. isotropy perturbation instability crossing declared tolerance,
3. reproducible adversarial sign-flip counterexample,
4. failed independent replay.

Then status must revert at least to `CLOSED_FULL_AXIOM_APPROVED_SB1` or lower,
with explicit incident packet.

---

## 5) For non-experts (PL, laika)

W praktyce: nie tylko „wymyśliliśmy” źródło selektora, ale też je
przetestowaliśmy pod obciążeniem.

Skutek: pełne domknięcie kernel-alone w wersji strict jest utrzymane,
ale z bezpiecznikiem — jeśli pojawi się mocny kontrprzykład, status wraca w dół.

To jest dokładnie styl uczciwej nauki: mocna teza + jawny warunek cofnięcia.
