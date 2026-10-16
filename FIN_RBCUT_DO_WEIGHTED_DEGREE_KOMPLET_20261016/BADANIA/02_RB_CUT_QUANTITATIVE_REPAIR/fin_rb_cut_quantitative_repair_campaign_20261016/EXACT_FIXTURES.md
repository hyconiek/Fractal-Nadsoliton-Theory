# Exact finite rational controls for the new cut/reconstruction theorem

Every finite example below is an **artificial graded reversible transport graph**, NOT the q12/A7 MT Gibbs count domain. Every number is checked with Python `fractions.Fraction`, by exhaustive subset cuts plus independent augmenting-path max-flow.

## 4 states, exact collective repair

Levels `0:(0)`, `1:(1,2)`, `2:(3)`; phases A=`{0,2,3}`, B=`{1}`. Allowed directed H edges `0->2, 2->3`, each reference conductance 1. Native-deleted phase-crossing edges F are `0->1, 1->3`, each conductance 1. Full native source is `s=(2,0,0,-2)`.

The sharp multiplier is `K_star=2`; flow `c_H(0->2)=2`, `c_H(2->3)=2` reconstructs **the exact full native source**, even though neither forbidden edge can be rerouted to its original endpoint by a positive H-path. K<2 fails a cut `S={0}`. The B phase has zero net source.

## 8 states, signed rather than additive correction is necessary

Levels: 0=`{0,1}`, 1=`{2,3}`, 2=`{4,5}`, 3=`{6,7}`. Phase A=`{1,2,4,6}`, B=`{0,3,5,7}`. Allowed native H edges:

- A: `1->2`, `2->4`, `4->6`, each reference conductance 2;
- B: `0->3`, `3->5`, `5->7`, each reference conductance 1.

Deleted phase-crossing F edges: `0->2` and `4->7`, each reference conductance 1. The exact native source `s=A_(H union F)c0` has a legal H-only solution: A flows `(2,1,2)` and B flows `(2,2,2)`. Its `K_star=2`.

But **the correction `delta=c_H-c0_H` decreases `2->4` by 1**. No purely nonnegative additive delta can exist: for `S={4,6}`, the deleted-edge divergence satisfies `sum_S b=+1`, while H has `C_H^+(S)=0`. It violates the additive-only cut criterion for every D. This explicitly prevents an erroneous negative verdict based on failed positive-only rerouting.

## 2 states, a genuine finite ideal obstruction (not FIN)

H has no edge; a native forbidden edge `0->1` of conductance1 yields `s=(1,-1)`. `S={0}` has positive source 1 and zero allowed outward capacity. Therefore `K_star=+infinity`, equivalently the complementary lower ideal has negative total source.

## Diverging relative cap can coexist with bounded absolute activity

In the 4-state example replace the two H native reference conductances 1 by `1/N`, leaving the two F conductances 1. Then the valid H flow has edge weights `1+1/N`; `K_star=N+1` diverges. With uniform vertex weights W=1 its total exit rates are <=2. This proves the reference-relative criterion is strictly sufficient, not necessary, for a uniform absolute kinetic cap, in the category of artificial graded networks.

## Audit limitations

Checks are algebraic/finite examples only, not the actual all-N q12 graph. No abstract toy discharges the H-ideal sign, uniform relative cut bound, uniform original vertex cap, or stochastic IBA limit on the MT domain.
