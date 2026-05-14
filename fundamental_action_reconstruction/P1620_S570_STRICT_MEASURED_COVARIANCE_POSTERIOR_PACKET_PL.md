# P1620 / S570 — Strict posterior with measured discriminant covariance

## Cel
Przejść z syntetycznego posterioru do posterioru opartego o macierz kowariancji
D1/D2/D3 (pomiarową lub referencyjną).

## Wejścia
- `generated/p1619_s569_strict_bayesian_posterior_sampler_summary.json`

## Wyjście
- `generated/p1620_s570_strict_measured_covariance_posterior_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Jawny model kowariancji i wpływu systematyk na posterior.
- Bez walidacji zewnętrznej.
