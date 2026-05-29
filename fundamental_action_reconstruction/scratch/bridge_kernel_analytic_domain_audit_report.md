# Scratch analytic-domain audit: legacy vs strict kernel

Status: analytic-domain tradeoff recorded; no bridge theorem.

- Legacy denominator is meromorphic with pole `z=-100`; strict eta is exact `9/5` with branch points `z=0, infinity` and finite cover `z=t^5`.
- `D_f-1=1.772588722239781` differs from strict `9/5` by `0.027411277760219` and is not a finite-sheet rational exponent.
- `9/5` is not the best tested rational approximation to `D_f-1`: all-tested-best flag `False`.
- Prior denominator-placement win replayed: `True`.
- No false pass: no kernel identity, no D_f-1→9/5 rationalization theorem, no physical-role transfer, no QW-2191 discharge, no ToE closure.
