# Scratch strict-alpha binary-rescale quotient probe

Status: quotient/canonicalization certificate for eta=9/5; no physical gauge theorem.

- Reduction rule: `(D,m,n)->(D/2,m,n-m)` preserves `2^n/D^m` when `D` is even and `n>=m`.
- Bounded scan: `56` rows through `k<=6`, `t<=8`; all reduce to odd representatives with steps equal to `k`.
- Canonical result: every `D=3*2^k,m=5*t,n=(8+5*k)*t` row reduces to `D=3,m=5*t,n=8*t`.
- Honest read: D=3 is canonical under binary-rescale quotienting, but the quotient/gauge principle itself is not derived.
- No false pass: no binary-rescale gauge theorem, no selector discharge, no QW-2191 discharge, no ToE closure.
