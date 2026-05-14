# P1621 / S571 — Publication-grade uncertainty budget

## Cel
Wyeksportować pełny budżet niepewności (stat+sys) dla posterioru strict
oraz plan reprodukowalnego notebooka.

## Wejścia
- `generated/p1620_s570_strict_measured_covariance_posterior_summary.json`

## Wyjście
- `generated/p1621_s571_publication_grade_uncertainty_budget_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Jawne rozbicie statystyczne/systematyczne per parametr kernela.
- Bez walidacji zewnętrznej.
- Status strict-core closure pozostaje `OPEN`, dopóki nie ma pełnych eksportów/witnessów/theoremów.

## Kontekst łańcucha teoretycznego
- Ten krok dotyczy wyłącznie warstwy niepewności na torze:
  `K_strict -> współczynniki -> pełny Lagrangian (L_SM + L_GR + sprzężenia) -> EOM`.
- Budżet niepewności nie zastępuje brakujących dowodów domknięcia strict-core.
