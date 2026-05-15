# P1718 S668 Strict Metric Tensor Backend Interface Contract Packet (PL)

Status: `P1718_EXECUTED_STRICT_TENSOR_BACKEND_INTERFACE_CONTRACT_NO_FALSE_PASS`  
As of: `2026-05-15`

## Cel

Po `P1717` zamienić blokadę rachunkową na konkretny kontrakt interfejsu backendu
liczącego residual metryczny.

## Co wyeksportowano

1. Dwie opcje backendu (componentwise_sympy / xAct-like).
2. Minimalny kontrakt I/O dla testu `EL_g - E_munu`.
3. Kontrakt eksportu obstrukcji przy wyniku niezerowym.
4. Decyzję wykonawczą: start od `componentwise_sympy`.

## Rygor

- strict-only,
- bez bridge do legacy,
- bez fałszywego pass,
- status globalny nadal `KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED`.

## Dla laika

To techniczna mapa narzędzia, które ma policzyć najtrudniejszy test grawitacyjny.
Zamiast ogólnego „brakuje rachunku”, mamy teraz dokładną specyfikację jak ten
rachunek wykonać i jak raportować wynik.

## Następny uczciwy krok (rekomendacja)

Zaimplementować runner `componentwise_sympy` i opublikować pierwszy wynik:
`residual-zero` albo jawny certyfikat obstrukcji.
