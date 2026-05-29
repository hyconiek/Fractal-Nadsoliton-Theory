# Scratch strict-alpha Fourier minimum-ripple convexity certificate probe

Status: integer convexity certificate for the phase/Fourier minimum-ripple selector; no physical selector theorem.

- Conditional theorem: for fixed `m,n`, Parseval reduces minimum non-DC Fourier ripple to minimizing `sum(e_i^2)`.
- Pairwise smoothing certificate: `(a,b)->(a-1,b+1)` for `a>=b+2` drops ripple by `2*N*(a-b-1)`.
- Target result: for `m=5, n=8`, the unique canonical minimizer is `[2, 2, 2, 1, 1]` and `q^5=256/243`.
- Bounded scan: `180` `(m,n)` cases and `5583` canonical ledgers checked; all had unique balanced minimizers: `True`.
- Honest read: this proves the integer convexity part only after the minimum-ripple premise is supplied; it does not derive the strict-core selector source.
- No false pass: no strict phase/Fourier selector theorem, no QW-2191 discharge, no ToE closure.
