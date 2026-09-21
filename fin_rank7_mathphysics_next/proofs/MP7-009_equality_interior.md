# MP7-009 — equality classification in the strictly positive interior

Assume `a3,a4,a5>0`.

For `b>0`, MP7-007 gives exact support

`L6=ker[f6]`,  `f6(m)=3m3+4m4+5m5 mod 6`.

The vectors

`(2,0,0)`, `(1,1,1)`, `(0,3,0)`

belong to `L6` and have determinant `-6`.  Since `f6` is surjective (`gcd(3,4,5,6)=1`), `L6` has index six in `Z^3`; therefore these vectors form a lattice basis.

For `b=0`, the support is

`L12=ker[f12]`, `f12(m)=3m3+4m4+5m5 mod 12`.

The vectors

`(4,0,0)`, `(1,1,1)`, `(0,3,0)`

have determinant `-12`; `f12` is surjective, so they form a basis of the index-twelve kernel.

Equality in MP7-008 means the phase character `chi_phi(m)=e^(im.phi)` is trivial on the relevant kernel.  Thus it factors through the cyclic image.  Every character of the image subgroup extends to the ambient cyclic group, giving

- for `b>0`: `phi_k = 2*pi*q*k/6 (mod 2*pi)`, `q in Z6`;
- for `b=0`: `phi_k = 2*pi*q*k/12 (mod 2*pi)`, `q in Z12`.

For positive `b`, these are exactly the even label translations (`l=2q`) preserving the alternating sign.  For `b=0`, they are arbitrary label translations.  Reflection sends `q` to `-q` and introduces no additional equality fields.

Thus in the strictly positive three-mode interior there are six equality fields at fixed `b>0` and twelve at `b=0`, before any accidental stabilizer from vanishing amplitudes (handled separately in MP7-010).
