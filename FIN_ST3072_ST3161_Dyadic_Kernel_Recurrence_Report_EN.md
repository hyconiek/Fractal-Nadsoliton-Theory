# FIN ST3072–ST3161 — Dyadic Kernel Recurrence

## Candidate

Use the literal strict formula on dyadic scales:

\[
q_n=K_{\rm strict}(2^n),
\qquad
K_{\rm strict}(d)=\frac{\cos(0.18575d+0.1625)}{1+d^{1.8}}.
\]

For \(\theta_0=0.3\), the resulting exact prism hierarchy converges to

\[
\boxed{\theta_\infty=0.32668864347458015.}
\]

The absolute cubic sum is

\[
\sum_n|\tanh q_n|^3=0.09106826887155.
\]

## Rigorous tail theorem

For \(d\ge1\),

\[
|K(d)|\le d^{-\eta},
\qquad
|\tanh K(d)|\le|K(d)|.
\]

Hence

\[
|\tanh K(2^n)|^3\le2^{-3\eta n}.
\]

The dyadic hierarchy is absolutely summable for every \(\eta>0\). Oscillatory
signs do not affect convergence. The extrapolated tail beginning at \(d=8\)
changes the declared limit by only about \(9.7\times10^{-8}\).

## Provenance obstruction

The frozen strict operator uses cyclic distances \(1,\ldots,6\). Values at
\(8,16,\ldots\) are evaluations of the analytic formula outside the exported
finite carrier. No current theorem identifies refinement layer \(n\) with
distance \(2^n\).

Thus the recurrence is stable but conditional:

\[
\text{formula evaluability}\neq\text{strict operator provenance}.
\]

## Scale-origin ambiguity

The generalized family

\[
q_n=K(a r^n),\qquad r>1,
\]

is summable for every positive \(a\) and every \(r>1\), but produces different
limits. For the dyadic ratio and \(\theta_0=0.3\), changing \(a\) from (0.5)
to (6) changes the limit from approximately (0.408) to (0.3000002).

Binary refinement makes (r=2) natural, but it does not prove that kernel
distance is multiplied by two at each informational layer. Neither (a) nor a
physical layer spacing is selected.

## Legacy comparison

Literal legacy dyadic sampling is also cubically summable, because its envelope
decays as (O(2^{-n})). Its early amplitudes and signs give substantially
different limits. This shared convergence class is not a legacy-to-strict
bridge and transfers no legacy physical role.

## Gate

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Final interpretation

The strict kernel shape can regularize a dyadic feedback hierarchy. The current
repository does not derive the infinite-distance extension, layer-to-distance
embedding, scale origin, refinement ratio or physical spacing. The highest
value next theorem is therefore a typed embedding from FIN refinement levels
to kernel distance.

No selector, physical layer scale, units, apparatus evidence, Standard Model,
gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
