# P1619 / S569 — Strict Bayesian posterior sampler (synthetic)

## Cel
Uruchomić syntetyczny sampler posterioru dla parametrów kernela strict
na podstawie discriminantów D1/D2/D3.

## Wejścia
- `generated/p1618_s568_strict_data_assimilation_protocol_summary.json`
- `generated/p1616_s566_strict_phenomenology_bounded_predictions_summary.json`

## Wyjście
- `generated/p1619_s569_strict_bayesian_posterior_sampler_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Jawna procedura posterioru i diagnostyki degeneracji.
- Bez walidacji zewnętrznej.
