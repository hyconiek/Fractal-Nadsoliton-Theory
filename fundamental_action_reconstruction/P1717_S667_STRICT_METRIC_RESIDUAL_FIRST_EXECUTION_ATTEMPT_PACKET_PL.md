# P1717 S667 Strict Metric Residual First Execution Attempt Packet (PL)

Status: `P1717_EXECUTED_STRICT_METRIC_RESIDUAL_FIRST_ATTEMPT_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po normalizacji konwencji (`P1716`) wykonać pierwszą realną próbę testu
`EL_g - E_munu` i jawnie wyeksportować wynik lub obstrukcję.

## Co wyeksportowano

1. Status próby wykonawczej: `BLOCKED_UNEXPANDED_TENSOR_CALCULUS`.
2. Konkretną listę członów blokujących rachunek.
3. Jawne rozdzielenie: co już gotowe vs czego brakuje.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To uczciwy raport z próby: wiemy już dokładnie, że problemem jest nie idea,
lecz trudny rachunek tensorowy wyższych członów krzywizny. Dzięki temu następny
krok może być precyzyjnie techniczny, a nie kolejne ogólniki.

## Następny uczciwy krok (rekomendacja)

Podłączyć backend do pełnej wariacji tensorowej i opublikować wynik:
`residual-zero` albo formalny certyfikat obstrukcji.
