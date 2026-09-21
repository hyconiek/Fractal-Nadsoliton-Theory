# MP7-041 — validated amplitude-to-phase continuation

Scientific state: **PROVED_INTERVAL_ASSISTED_EXACT_SYMMETRY_CONTINUATION**.

## Selected branch

Use full-phase root ID 14 from the accepted R7P-092 fixed-amplitude catalog.  Its phase lock is

`(phi3,phi4,phi5) = (pi/2, 2 pi/3, 5 pi/6)`

with the fixed-fixture sign convention `r3,r4,r5>0`, `z6<0`.

This is not an arbitrary numerical root: translation `j -> j+1` changes the aligned positive-`z6` field into exactly this phase lock while flipping the alternating harmonic.  Consequently, for **every** amplitude tuple with the same signs, the phase gradient of `log Z` vanishes at this lock.  Hence the local continuation map is stronger than a merely implicit map:

`phi(a) = (pi/2,2pi/3,5pi/6)` exactly.

## Certified amplitude box

Around the exact-decimal fixture

`a0=(0.1131879146, 0.1698528641, 0.2269339093, -0.3380663037)`

use the independent relative box of radius `0.001` (0.1%).  Its bounds are

- `r3 in [0.1130747266854, 0.1133011025146]`,
- `r4 in [0.1696830112359, 0.1700227169641]`,
- `r5 in [0.2267069753907, 0.2271608432093]`,
- `z6 in [-0.3384043700037,-0.3377282373963]`.

All signs are preserved.

Direct interval evaluation of the full log-partition phase Hessian on this box, followed by LDL on the **Phi** angular block `-H_phi_phi(log Z)`, gives positive pivots

1. `[5.93748118897e-4, 7.18505111299e-4]`,
2. `[4.41932102218e-4, 7.09184815707e-4]`,
3. `[3.61520112818e-5, 4.57883416279e-4]`.

Thus the phase block is uniformly nonsingular and positive for `Phi`; the IFT condition is paid throughout the box.  The exact symmetry identity already supplies the continuation itself.

## Cross block and Schur complement

Because `grad_phi log Z` is identically zero as a function of the amplitudes on this symmetry-locked family,

`H_{a phi}(log Z)=0`

exactly.  Therefore the phase-eliminated amplitude Hessian does not need a numerical inverse/cancellation:

`H_eff = H_aa(Phi)`.

In the complex-amplitude convention of the phase fixture,

`||theta||^2 = 2 r3^2/lambda3 + 2 r4^2/lambda4 + 2 r5^2/lambda5 + z6^2/lambda6`,

so at the center

`H_eff(g)=Q/g-K_aa`,

where

`Q=diag(2/lambda3,2/lambda4,2/lambda5,1/lambda6)`

and `K_aa` is the covariance Hessian of `log Z` with respect to `(r3,r4,r5,z6)`.

The serialized result contains both matrices and sample spectra for `g=3.5,3.7,4,5`.

## Connection to MP7-040

At the frozen fixture, the four componentwise radial gain candidates are numerically

`5.0624496338, 5.0624496875, 5.0624496750, 5.0624496775`.

Their closeness explains why the symmetry-locked phase root looks almost like a full equilibrium in ordinary floating arithmetic.  It does **not** override MP7-040: its outward interval check proves that the required radial ratios are not exactly equal at the frozen fixture.

The new conclusion is different and useful: nearby full equilibria on this symmetry family can be searched entirely in the four radial amplitudes, because the phase coordinates are fixed exactly by symmetry.  No frozen-phase census needs to be extrapolated by assumption.

## Nonconclusions

This continuation certificate does not prove that a nearby radial/full equilibrium exists, does not classify all amplitude continuations of the other 59 phase roots, and does not imply local or global minimality in the amplitude directions.
