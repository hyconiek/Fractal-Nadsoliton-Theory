# MULTITIME-SEMIGROUP-FINGERPRINT-298
## Multi-time rate drift is an exact memory fingerprint, while the accepted N=6 FIN localized process converges to a single late-time 12-state generator after a short microscopic boundary layer

Date: 2026-09-27

Status:
- exact theorem for finite reversible microscopic dynamics projected onto one symmetry-adapted scalar mode;
- secondary numerical analysis of the accepted exact N=6 microscopic correlation artifact from report 229;
- exact algebra for the D12 shell reconstruction;
- finite-sample Chernoff bounds under the same illustrative symmetric readout-noise model as report 297;
- no new full microscopic eigensolve was performed in this report, so the N=6 source correlations retain the epistemic status of their accepted upstream artifact.

Source Git blob SHAs used for provenance:

- `micro_to_12_semigroup_error_N6.json`: `fd9425f16a0aa20e293a231ea5f6a171dcac8c7e`
- `localized12_MZ_fourier_N6.json`: `0fb7a9e9daad2ea1e2f2cc4d58dc2c247ad297cb`
- `localized12_exact_mode_validation_N6.json`: `d82d9cd827a39b952c45c0837454728d60c17094`

## 1. Question

Report 297 proved that one localized preparation plus the seven-bin scalar readout

\[
Y_1(J)=\cos(2\pi J/12)
\]

contains all one-time information available in the declared reflection-symmetric D12 effective class. Therefore adding another one-time observable cannot test whether the 12-state law is genuinely Markovian.

The new question is genuinely multi-time:

> after fitting the resolved dynamics at one time or in one time window, does the same generator predict the process at other times?

A naive test such as `log C_k(t)/t = constant` is insufficient because report 229 already found a fixed initial-slip residue

\[
C_k(t)\simeq Z_k e^{\lambda_k t},\qquad Z_k<1,
\]

after the microscopic memory boundary layer. The factor `Z_k` would spuriously make `log C/t` time dependent even when the remaining late dynamics is a perfect single exponential.

The correct diagnostic must eliminate any fixed multiplicative residue.

## 2. Exact reversible spectral theorem

Let `S` be the self-adjoint reversible microscopic generator in the equilibrium Hilbert space, and let `b_k` be a normalized resolved Fourier observable. Since `-S` is positive semidefinite, the spectral theorem gives

\[
\boxed{
C_k(t)=\langle b_k,e^{tS}b_k\rangle
      =\sum_\alpha w_{k\alpha}e^{-r_{k\alpha}t}
}
\]

with

\[
w_{k\alpha}\ge0,\qquad \sum_\alpha w_{k\alpha}=1,\qquad r_{k\alpha}\ge0.
\]

Hence `C_k` is a positive mixture of exponentials.

Define the instantaneous effective relaxation rate

\[
r_k(t)=-\partial_t\log C_k(t).
\]

With the tilted spectral weights

\[
\widetilde w_{k\alpha}(t)
 =\frac{w_{k\alpha}e^{-r_{k\alpha}t}}{C_k(t)},
\]

we get exactly

\[
\boxed{
r_k(t)=\mathbb E_t[r]}
\]

and

\[
\boxed{
r'_k(t)=-\operatorname{Var}_t(r)\le0.}
\]

Equivalently,

\[
\boxed{\frac{d^2}{dt^2}\log C_k(t)=\operatorname{Var}_t(r)\ge0.}
\]

So `log C_k(t)` is convex.

This supplies an exact interpretation:

- a single exponential has zero spectral variance and a constant rate;
- any nonzero rate drift witnesses more than one microscopic decay rate in the resolved observable;
- the rate drift is therefore a direct memory / unresolved-spectrum fingerprint.

The theorem does not require fitting a Mori-Zwanzig kernel.

## 3. Three-time residue-free semigroup statistic

For `a<b<c`, define the two divided slopes

\[
s_k(a,b)=\frac{\log C_k(b)-\log C_k(a)}{b-a},
\]

\[
s_k(b,c)=\frac{\log C_k(c)-\log C_k(b)}{c-b},
\]

and the drift

\[
\boxed{
\Delta_k(a,b,c)=s_k(b,c)-s_k(a,b).
}
\]

Convexity proves

\[
\boxed{\Delta_k(a,b,c)\ge0.}
\]

For any pure exponential with arbitrary fixed residue,

\[
C(t)=Z e^{\lambda t},
\]

both divided slopes equal `lambda`, so

\[
\boxed{\Delta=0.}
\]

Thus this statistic removes the fixed initial-slip amplitude `Z`. A positive value is not produced merely by `Z<1`; it requires continued multi-rate structure across the selected time windows.

## 4. Readout-noise cancellation

Keep the symmetric readout channel from report 297:

\[
P(J_{\rm rep}=j|J=i)=
\begin{cases}
1-\eta,&j=i,\\
\eta/11,&j\ne i.
\end{cases}
\]

For every nonzero Fourier mode,

\[
\widetilde C_k(t)=a_\eta C_k(t),
\qquad
 a_\eta=1-12\eta/11.
\]

Therefore

\[
\log\widetilde C_k(t_b)-\log\widetilde C_k(t_a)
=
\log C_k(t_b)-\log C_k(t_a),
\]

and the factor `a_eta` cancels exactly.

Hence, in the infinite-sample limit,

\[
\boxed{
\Delta_k\text{ is invariant under any time-independent symmetric readout error }\eta.
}
\]

This is stronger than the one-time test of report 297, where the noise level reduces the Chernoff exponent. Here the mean rate-drift itself does not depend on `eta`; noise only enlarges finite-sample uncertainty.

## 5. Exact N=6 microscopic replay: mode k=2 is the strongest early memory carrier

Using the accepted exact projected correlations from report 229, the interval slopes for `k=2` are:

| interval | divided slope |
|---|---:|
| 0.10 -> 0.25 | -0.0966295872 |
| 0.25 -> 0.50 | -0.0674347946 |
| 0.50 -> 1 | -0.0489278676 |
| 1 -> 2 | -0.0398563373 |
| 2 -> 4 | -0.0360868863 |
| 4 -> 8 | -0.0350595733 |
| 8 -> 16 | -0.0349625033 |
| 16 -> 32 | -0.0349610440 |
| 32 -> 64 | -0.034961043177 |

They increase monotonically, exactly as the positive-spectral-measure theorem predicts.

The exact slow microscopic eigenvalue in this symmetry sector is

\[
\lambda_{2,\rm exact}
=-0.0349610431769811.
\]

The late 32 -> 64 slope differs from it by only

\[
3.5\times10^{-14}
\]

relatively for `k=2`.

Across all six nontrivial modes, the largest relative difference between the 32 -> 64 incremental slope and the corresponding exact microscopic slow eigenvalue is

\[
\boxed{2.61\times10^{-12}.}
\]

This is a much stronger statement than saying the late curves merely look exponential.

## 6. The memory curvature collapses rapidly

The largest three-time drift over the six modes is:

| triple `(a,b,c)` | max Delta | strongest mode |
|---|---:|---:|
| (0.10, 0.25, 0.50) | 2.91948e-2 | 2 |
| (0.25, 0.50, 1) | 1.85069e-2 | 2 |
| (0.50, 1, 2) | 9.07153e-3 | 2 |
| (1, 2, 4) | 3.76945e-3 | 2 |
| (2, 4, 8) | 1.02731e-3 | 2 |
| (4, 8, 16) | 9.70699e-5 | 2 |
| (8, 16, 32) | 1.95190e-6 | 5 |
| (16, 32, 64) | 4.53157e-9 | 5 |

Thus the projected process is decisively not a single Markov exponential near `t=0`, but its multi-rate curvature becomes negligible on the late window.

This resolves an apparent tension between reports 218/229 and a strict semigroup test:

\[
\boxed{
\text{microscopic projected process}
=\text{short memory boundary layer}
+\text{very accurate late single-pole dynamics}.
}
\]

The 12-state generator should therefore be interpreted as a controlled late-time effective generator, not as an exact microscopic semigroup starting at `t=0`.

## 7. The late incremental spectrum reconstructs the exact 12-state shell generator

Use the six incremental rates from the window 32 -> 64 and invert the exact shell matrix

\[
\lambda_k
=2\sum_{d=1}^{5}q_d[\cos(2\pi kd/12)-1]
+q_6[(-1)^k-1].
\]

As in report 297, this matrix has

\[
\det M=-3456\sqrt3\ne0.
\]

The shell rates reconstructed from the late incremental slopes agree with the independently stored exact spectral shell rates of report 214 with maximum relative difference

\[
\boxed{3.72\times10^{-11}.}
\]

Therefore the multi-time test does more than reject exact short-time Markov closure: it independently identifies the time window in which the same six-rate D12 generator becomes dynamically exact to numerical precision for the resolved slow poles.

## 8. Relation to the optimal one-time fingerprint of report 297

For N=6,

\[
\rho=0.02261891215447.
\]

The report-297 optimum

\[
\tau_*=\rho t_*=0.5427059873
\]

corresponds to

\[
\boxed{t_*\approx23.9935}
\]

in microscopic refresh-time units.

This lies deep in the regime where the rate curvature is already tiny: the 8 -> 16 -> 32 drift is only of order `10^-6`.

Therefore the two fingerprints probe different physics:

- **297:** resolved shell structure around `tau ~ 0.54`, after short memory has mostly died;
- **298:** early-time non-semigroup / memory structure around `tau << 1`.

A single measurement time cannot optimize both goals.

This is useful experimentally because it prevents a false expectation that the strongest FIN-vs-comparator shell fingerprint should also be the strongest memory test.

## 9. Same physical readout is enough

No new observable is required.

The seven-bin `Y_1` histogram from report 297 reconstructs the reflection-symmetric 12-label distribution and therefore every cosine Fourier coefficient, including the strongly memory-sensitive `k=2` mode.

Thus one apparatus can perform both tests by changing only the evolution time:

1. prepare localized `J=0`;
2. read the same scalar `Y_1`;
3. use several early and late time bins;
4. reconstruct Fourier coefficients;
5. test slope constancy and the resolved shell fingerprint separately.

## 10. Concrete two-window validation protocol

A particularly clean N=6 protocol is:

1. estimate the late incremental generator from the ratio of Fourier coefficients between `t=32` and `t=64`;
2. this automatically cancels the fixed pole residue `Z_k`;
3. invert the six incremental slopes to `q_d`;
4. without refitting, extrapolate that late Markov generator back to `t=1`;
5. compare the predicted and observed seven-bin `Y_1` histogram.

At 5% symmetric readout noise, the exact microscopic projected distribution at `t=1` and the late-generator prediction have

\[
C_{\rm Chernoff}=9.1593954\times10^{-4}.
\]

Conditioned on an independently calibrated late generator, the standard equal-prior Chernoff bound gives a sufficient validation-shot count

\[
\boxed{M=2514}
\]

for

\[
P_e\le5\%,
\]

and

\[
M=4272
\]

for the corresponding 1% bound.

Noise robustness for the same frozen early/late protocol is:

| eta | Chernoff C | shots for <=5% bound |
|---:|---:|---:|
| 0% | 0.00219309 | 1050 |
| 5% | 0.000915940 | 2514 |
| 10% | 0.000568821 | 4048 |
| 20% | 0.000292839 | 7863 |

These counts do **not** include uncertainty in estimating the late generator. They are validation counts after an independent late-time calibration.

## 11. Direct rate-drift test without a fitted comparator

A model-independent alternative uses only mode `k=2` at three times `(0,1,32)`.

The measured slopes are

\[
s_2(0,1)=-0.0692851038,
\]

\[
s_2(1,32)=-0.0352046818,
\]

so

\[
\boxed{\Delta_2=0.0340804221.}
\]

Under 5% symmetric readout noise, a delta-method normal approximation with one-sided `alpha=0.05` and 80% power gives approximately:

- equal allocation: 1771 shots at each of the three times, about 5312 total;
- variance-optimal allocation: about 4034 total shots, distributed approximately `34.1% / 55.9% / 10.0%` across `t=0 / 1 / 32`.

This is an asymptotic normal power estimate, **not** an exact finite-sample guarantee. It is included to show the scale of the experiment and to demonstrate that the rate-drift statistic remains operational under readout noise.

## 12. What 298 establishes

Inside the accepted N=6 microscopic lane:

- reversible projected correlations have a positive spectral representation;
- their logarithmic curvature is an exact unresolved-spectrum/memory diagnostic;
- the `k=2` mode carries the strongest early curvature on the stored grid;
- the effective rate moves monotonically toward the exact microscopic slow eigenvalue;
- late 32 -> 64 slopes reproduce all six exact slow eigenvalues to relative error below `2.7e-12`;
- inverting those slopes reproduces the exact six shell rates to below `3.8e-11` relative error;
- the one-time `Y_1` sensor from 297 already contains all data required for the multi-time test;
- fixed symmetric readout attenuation cancels exactly from divided logarithmic slopes;
- the shell fingerprint and the memory fingerprint occupy parametrically different time windows.

## 13. What 298 does not establish

- N=6 does not prove an N-uniform memory theorem;
- the collapse of curvature at N=6 does not prove that the boundary layer stays O(1) as N -> infinity;
- the localized 12-state labels are not thereby physical spatial sites;
- `rho` is still an internal time calibration, not an SI clock;
- no apparatus or measured noise level is supplied;
- finite-sample counts assume independent preparations;
- the Chernoff counts condition on an independently calibrated late generator;
- the result does not distinguish FIN from every possible non-Markov model, only tests the declared semigroup property of this reduced lane.

## 14. New research implication

The strongest next question is no longer whether residual memory exists at N=6; it does, and its decay into the late generator is now quantitatively resolved.

The next proof-grade target should be:

### 299 — MEMORY-WINDOW-SCALING

Determine whether the dimensionless width of the microscopic boundary layer shrinks relative to the metastable clock as N increases.

A useful target quantity is

\[
\epsilon_{\rm mem}(N)
=\rho_N\,t_{\rm settle}(N;\varepsilon),
\]

where `t_settle` is defined by a residue-free multi-time criterion such as

\[
\max_k\Delta_k(t,2t,4t)\le\varepsilon |\lambda_{k,\rm slow}|.
\]

If

\[
\epsilon_{\rm mem}(N)\to0,
\]

then the memory-aware 12-state generator gains a genuine scale-separation theorem rather than only finite-N validation.

If it does not tend to zero, the memory kernel must remain explicit in the large-N effective theory.

## Verdict

P1-298 is positive and sharpens the effective-theory interpretation.

The microscopic-to-12-state reduction is **not** an exact Markov semigroup from time zero. The exact reversible spectral representation forces a measurable positive rate drift during the initial memory layer.

But that drift collapses rapidly, and the late 32 -> 64 incremental spectrum at N=6 reconstructs the exact slow six-shell generator essentially to numerical precision.

The strongest current statement is therefore

\[
\boxed{
\text{exact microscopic FIN}
\to
\text{short multi-rate memory layer}
\to
\text{single late D12 Markov generator}
}
\]

rather than either of the two overstatements

\[
\text{exact microscopic FIN}=\text{12-state Markov chain from }t=0
\]

or

\[
\text{memory prevents a controlled 12-state Markov description}.
\]
