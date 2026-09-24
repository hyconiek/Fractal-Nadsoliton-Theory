# DYN-002 derivation note

For one coordinate, write `Delta_+=E(x+eps)-E(x)` and `Delta_-=E(x-eps)-E(x)`. With

`q_+=M0/(beta eps^2) exp(-beta Delta_+/2)`,
`q_-=M0/(beta eps^2) exp(-beta Delta_-/2)`,

expand

`Delta_+=eps E'+eps^2 E''/2+O(eps^3)`,
`Delta_-=-eps E'+eps^2 E''/2+O(eps^3)`

and `f(x+/-eps)-f(x)=+/-eps f'+eps^2 f''/2+O(eps^3)`. The `1/eps` terms cancel and

`q_+[f(x+eps)-f(x)] + q_-[f(x-eps)-f(x)]`
`= -M0 E' f' + (M0/beta) f'' + O(eps)`.

The multi-coordinate formula follows by summation.

Detailed balance is exact before taking the limit. Consequently the stationary density is proportional to `exp(-beta E)` whenever normalizable and boundary conditions are compatible.

The PHA-001 gap controls only transverse phase locking on its declared nonzero-amplitude box. It does not provide the missing radial spectral gap or eliminate the REF-002 response memory, so a full Markov slow-manifold theorem is not licensed here.
