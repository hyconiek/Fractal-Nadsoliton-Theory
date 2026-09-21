# MP7-038 — scoped C4-to-X7 stationary index bridge

At an aligned field (`phi=0`), the Gibbs distribution is reflection-even under `j -> -j`.  The four C4 features (three cosines plus the alternating cosine) are reflection-even, while the three sine features are reflection-odd.  Their cross covariances vanish exactly.  The Euclidean quadratic term is also block diagonal. Therefore the full Cartesian mediator Hessian splits

`H7 = H_even4 direct_sum H_odd3`.

Now assume:

1. the point is a full stationary point;
2. `s3,s4,s5>0` (interior in the three paired amplitudes);
3. after a label translation it is aligned with nonnegative alternating amplitude;
4. `0<g<=250/67`.

MP7-037 gives `H_odd3>0`.  MP7-005/006 gives `index_negative(H_even4)<=1` from Target P. Hence

`index_negative(H7) <= 1`

for this **interior aligned stationary family** throughout the Target-P gain window.

At `g=250/67`, Target P is non-strict, so the even block may have a zero mode. The odd block remains strictly positive under the interior premise.

## Why this does not resurrect the refuted unrestricted theorem

- The accepted nonstationary index-two witness is not stationary, so the polar gradient term need not vanish there.
- The accepted exact `g=5` stationary boundary counterexample lies outside `g<=250/67` and has missing amplitudes, so MP7-037's interior hypothesis does not apply.
- MP7-011 aligns **global minimizers**, not arbitrary stationary points. Thus this theorem does not exhaust nonaligned stationary families.

Boundary stationary roots require separate support/transverse analysis. MP7-013/014 prove that in the Target-P gain window no nonzero **aligned** boundary stationary support exists; this strengthens the aligned branch statement but still does not prove every stationary point is aligned.
