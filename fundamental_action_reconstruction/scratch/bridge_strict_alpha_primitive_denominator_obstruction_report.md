# Scratch strict-alpha primitive denominator obstruction probe

Status: primitive-denominator certificate for eta=9/5; no denominator theorem.

- Wide model: `q_model=(2^n/D^m)^(1/m)`; exact target requires `D=3*2^(n/m-8/5)`.
- Exact integer family: `D=3*2^k`, `m=5*t`, `n=(8+5*k)*t`; scanned `k<= 6`, `t<= 8` with `56` exact rows.
- Primitive branch: `k=0` recovers `D=3,m=5*t,n=8*t`; nonprimitive branches `k>0` are hidden-binary-rescale exact families.
- Honest read: denominator 3 is primitive inside the model, but the primitive/no-rescale convention is still not derived.
- No false pass: no denominator theorem, no selector discharge, no QW-2191 discharge, no ToE closure.
