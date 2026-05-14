# P1618 / S568 — Strict data-assimilation protocol

## Cel
Zmapować discriminanty (D1/D2/D3) na posteriory parametrów kernela strict
(omega, phi, beta, eta) z jawnymi założeniami o błędach.

## Wejścia
- `generated/p1617_s567_strict_experimental_discriminants_summary.json`
- `generated/p1616_s566_strict_phenomenology_bounded_predictions_summary.json`

## Wyjście
- `generated/p1618_s568_strict_data_assimilation_protocol_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Jawny model błędu i identyfikowalności dla każdego parametru.
- Bez walidacji zewnętrznej.
