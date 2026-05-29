# Scratch strict-alpha selector rational threshold certificate probe

Status: exact rational selector certificate for eta=9/5; no strict selector theorem.

- Rational rule: for `gamma=p/q`, balanced wins iff `3^q > 2^(p+q)`; rational tie is impossible by prime factorization.
- Bounded scan: `1261` reduced rationals with denominator <= `64`; safe `736`, fail `525`, tie `0`.
- Best safe in scan: `31/53`; first fail in scan: `24/41`; `gamma_c=0.584962500721`.
- Common certificates: `1/2` is safe; `3/5` is already above threshold and flips to `(3,2,1,1,1)`.
- No false pass: no strict derivation of a safe rational gamma, no QW-2191 discharge, no ToE closure.
