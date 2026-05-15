# P1764 / S714 — Jawny eksport kowariantny nonproxy: `E_A^μ`, `E_H`

Status: `P1764_S714_EXPLICIT_COVARIANT_NONPROXY_OPERATOR_EXPORT_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Podnieść stan po `P1763` (blokada template-only) do jawnego poziomu operatorowego:

```text
K_strict -> coefficient map -> full non-skeleton L_total -> explicit covariant E_A^μ, E_H
```

bez legacy bridge i bez proxy-skrótów.

## Co zostało dowiedzione (technical progress)

Wyeksportowano jawne formy kowariantne nonproxy (operator-level):

1. `E_A^μ := (1/sqrt(-g)) δS/δA_μ`:
   - `∇_ν(Z_a F_a^{νμ})`
   - `+ 2 χ_RG ∇_ν(R F_a^{νμ})`
   - `+ J_matter,a^μ + J_CT,a^μ = 0`
2. `E_H := (1/sqrt(-g)) δS/δH^†`:
   - `D_μD^μH + μ_H^2H + 2λ_H(H^†H)H + ξ_HR R H + λ_{φH} φ^2 H`
   - `+ d(CT_1loop)/dH^†`
   - `- Yukawa sources = 0`

To jest jawny postęp względem template-only: formuły są już fizycznie zakotwiczone
w `L_gauge/L_higgs/L_yukawa/L_scalar_phi/L_mix`.

## Co nadal jest OPEN

1. Brak pełnej postaci indeksowo-komponentowej dla wspólnej rodziny teł 4D.
2. Brak uruchomionego testu H1 z różnicą:
   `δE_A^μ/δH - δE_H/δA_μ`.
3. Brak bramkowania theorem-level dla:
   - Bianchi/Ward consistency,
   - BRST nilpotency,
   - Cutkosky unitarity,
   - renormalization/counterterm closure,
   - background-independence global consistency.

## Ryzyka false-pass

1. Mylenie „jawny operator kowariantny” z „pełny export komponentowy”.
2. Promocja weak-form readiness do strict-local PASS bez residual witness.
3. Ukrywanie wkładu `CT_1loop` bez jawnego śladu wariacyjnego.

W tym pakiecie **nie** wydano żadnego `PASS_ZERO`.

## Następny uczciwy krok (highest-value bottleneck)

1. Rozwinąć oba eksporty (`E_A^μ`, `E_H`) do pełnej postaci komponentowej
   na tej samej rodzinie teł i tej samej konwencji indeksowej.
2. Uruchomić test H1 4D:
   `δE_A^μ/δH - δE_H/δA_μ` z werdyktem tylko `PASS_ZERO` albo `OBSTRUCTION`.
3. Po H1, przejść na kanał metryczny (`EL_g - E_{μν}`) bez zmiany klasy teł.

## Dla laika

To jak przejście z projektu technicznego do pełnych równań do policzenia:
nie mamy już samych „szablonów”, tylko konkretne równania fizyczne. Ale żeby
uczciwie potwierdzić teorię, trzeba je jeszcze rozpisać do postaci, którą da się
bezpośrednio sprawdzić liczbowo/symbolicznie w każdym składniku.
