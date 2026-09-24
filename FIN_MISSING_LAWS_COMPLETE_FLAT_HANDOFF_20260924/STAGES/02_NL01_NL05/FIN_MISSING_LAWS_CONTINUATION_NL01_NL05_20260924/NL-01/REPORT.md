# NL-01 — multipole cutoff source

## Result

`NO_GO_L2_NOT_SELECTED_BY_CURRENT_COVARIANCE_COMPOSITION_CLASS`.

For a q-point circular carrier and integer L with q>2L, let `C_{q,L}` contain
monopole plus the real sine/cosine harmonics k=1,...,L.  DFT orthogonality gives

`rank C_{q,L} = 2L+1`,

so conditioning a full covariance that is nondegenerate outside its constant
mode gives an effective rank

`q-(2L+1)`.

L=2 therefore gives rank 7 at q=12, but L=0,1,3,... define equally covariant
families.  Rotation and reflection preserve every one of these row spaces.
Applying the same declared neutrality law independently to composed components
also preserves every L.  Hence the current symmetry and composition principles
select the *hierarchy*, not its cutoff.

The rule is also not rescued by ordinary locality: projecting out a finite set
of global Fourier harmonics is nonlocal for every L>=1.  Locality therefore does
not distinguish L=2 from neighboring choices.

## Consequence

The mediator proposal remains a useful conditional source candidate because it
has held-out predictions (`rank=q-5`, and conditional covariance rather than a
top-seven spectral cut).  But the rank-seven question is now reduced to the
new source question: **what independently selects multipole order L=2?**

A successful future source must predict the cutoff before seeing q=12/rank 7.
