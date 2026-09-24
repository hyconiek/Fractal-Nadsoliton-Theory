# LAW-01 — geometric neutrality as a rank-source candidate

**Disposition:** `PASS_AS_CONDITIONAL_SOURCE_CANDIDATE__NOT_STRICT_SOURCE`.

Let the full strict circulant Laplacian `A` be the covariance on the zero-sum mediator sector (equivalently, `A=BWB^T` can arise as the covariance of the divergence of independent Gaussian edge flows). Impose only dipole and traceless-quadrupole neutrality, represented by the real `k=1` and `k=2` Fourier rows `C`. Gaussian conditioning gives

`A_eff = A - A C^T (C A C^T)^(-1) C A`.

For the current q=12 strict operator this has rank **7** and agrees with the existing top-seven `A7` truncation to Frobenius residual `9.594e-15`. This is exact structurally for a circulant operator because the constraint row-space is precisely the invariant `k=1,2` eigenspace.

The important new content is outside q=12: the rule predicts `rank=q-5`, giving 13 at q=18 and 19 at q=24. Under a 20% perturbation of one edge, the conditioned covariance and a top-seven spectral truncation differ by relative Frobenius norm **4.0347%**; the conditioned law still satisfies the moment constraints (`1.44e-15`) while the top-seven truncation does not (`0.409`). Thus the proposal is falsifiably different from “keep the seven largest eigenvalues”.

## Source status

This does **not** derive the full strict `A`, Gaussianity, or physical mediator semantics. More importantly, neutrality through multipole order `L` gives a whole family `rank=q-(2L+1)`; selecting `L=2` is the remaining source datum. The proposal therefore compresses the old “why rank seven?” problem into the sharper “why exact neutrality through quadrupole and not higher/lower?” problem. It is a viable conditional parent law, not a strict FIN derivation.
