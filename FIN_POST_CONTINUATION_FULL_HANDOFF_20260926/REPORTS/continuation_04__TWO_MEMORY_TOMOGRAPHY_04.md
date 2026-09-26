# TWO-MEMORY-TOMOGRAPHY-04 — first-order full hidden-sector tomography

Status: **FIRST-ORDER CONDITIONAL IDENTIFIABILITY THEOREM**.

Arbitrary weak coupling can inject the tree state into the four hidden simplex
directions k=1,2.  A single k=5 quartic is insufficient because a hidden k=2
coupling can imitate its old nonzero/null-control pattern.

Use instead the preregistered mixed retained probe

f_t = A c_3 + t B c_4.

For a hidden cosine direction, the exact uniform-state quartic jump functional
has coefficients

K_1(t)=3 A B t (10 A^2+9 B^2 t^2),
K_2(t)=9 A^2 B^2 t^2.

Thus the odd-in-t part isolates k=1 and the even-in-t part isolates k=2.
Translate the probe through all twelve label shifts.  For each l=1,2 the
12-by-2 design matrix [cos(2 pi l a/12), sin(2 pi l a/12)] has Gram matrix 6 I
and rank 2.  Discrete Fourier projection therefore reconstructs independently
the cosine and sine amplitudes of both hidden sectors.

An exact symbolic fixture with hidden coefficients (a1,b1,a2,b2)=(2,3,5,7)
and t=1/3 gives K1=11, K2=1 and recovers:

odd k=1 projection = (22,33)=11(2,3),
even k=2 projection = (5,7)=1(5,7),

with exactly zero leakage into the wrong Fourier sector.  Uniform averaging
over the 12 shifts kills every first-order hidden contribution.

Because the intrinsic declared heat-bath dynamics is D12-equivariant, its
translation-invariant contribution lies in the zero shift harmonic and does
not contaminate the nonzero k=1,2 tomography harmonics.  The original k=5
phase scan can therefore be retained for the intrinsic constant quartic
coefficient while the mixed-probe scan diagnoses external hidden coupling.

Boundary: this is first order in the declared weak hidden-state coupling.
Angle-independent O(epsilon^2) terms or additional mechanisms in the same
zero harmonic can still confound the intrinsic constant component.
