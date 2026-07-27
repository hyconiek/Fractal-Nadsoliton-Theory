# Release 10.13 — FIN Programs 125–137

## Trace Selection, Natural Localizers, Fractional Physics, and Operational Sources

### Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 1.0.0
- **Publication date:** 2026-07-27
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Related repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Release 10.13 executes thirteen mathematical and operational research
programs. It tests whether the homological–character fibre built in Release
10.12 actually forces the strict damping exponent, strengthens the fractional
limit, and constructs the smallest missing measurement and calibration
objects.

The five fibre dimensions are

\[
(1,2,2,2,2).
\]

Even after granting the enlarged abstract symmetry
\(\operatorname{Aut}(U(12))\cong S_3\), the sector orbits are

\[
\{2\},\qquad\{3\},\qquad\{5,7,11\}.
\]

Every trace invariant under this enlarged symmetry therefore has the form

\[
(x,y,z/3,z/3,z/3),\qquad x+y+z=1,
\]

and the exponent readout is

\[
\eta=2-x.
\]

The uniform trace gives \(9/5\), but even this enlarged symmetry does not force
it. Changes of generator of \(C_{12}\) fix each multiplication map \(m_p\)
individually and leave an even larger trace simplex. The normalized Hilbert
trace gives \(17/9\). This is a new trace-selection no-go theorem.

The fibre assignment is nevertheless upgraded to a natural,
presentation-independent localizer on \(C_{12}\). A continuous
finite-frequency enclosure supports the \(4/5\)-fractional symbol law on
\([10^{-3},2\cdot10^{-2}]\), with maximum conservative relative-remainder
upper bound \(2.2647\%\). The result uses analytic tail and frequency-cell
bounds but guarded ordinary floating-point FFT rather than formal interval
arithmetic.

The same fractional generator supports unitary wave evolution in \(L^2\) and
stable diffusion. Pointwise wave records, however, require a UV resolution:
the standard dyadic curvature estimate grows as
\(t^{-1/2}\Lambda^{3/5}\) and is not ultraviolet summable.

The minimal physical calibration

\[
x_{\rm phys}=\ell x,\qquad
t_{\rm phys}=\tau t,\qquad
H_{\rm phys}=(\hbar/\tau)A
\]

has rank three. Length, time, and action remain independent external
calibrations. A nonparametric finite-memory apparatus estimator, an exact
projective crossover flow, an amplitude quotient, a trace-parameterized
damping completion, and a signed state receiver are constructed. The receiver
detects opposite prepared branches but does not choose one, so QW-2191 remains
open.

## Executed programs

1. **Program 125:** invariant trace classification and \(9/5\) no-go theorem.
2. **Program 126:** natural presentation-independent \(C_{12}\) fibre
   localizer.
3. **Program 127:** continuous finite-window fractional-symbol enclosure.
4. **Program 128:** quantitative finite-step stable comparison.
5. **Program 129:** fractional-wave UV obstruction and cutoff family.
6. **Program 130:** rank-three dimensional calibration object.
7. **Program 131:** nonparametric apparatus process tomography.
8. **Program 132:** exact projective local/fractional crossover flow.
9. **Program 133:** bounded phase/frequency source test and rejection.
10. **Program 134:** amplitude projectivization and role-loss theorem.
11. **Program 135:** trace-parameterized conditional damping bridge.
12. **Program 136:** signed operational state receiver.
13. **Program 137:** external-data intake audit.

## Principal scientific results

- The localizer is natural, but its state is not selected.
- Traces invariant under enlarged \(\operatorname{Aut}(U(12))\) form a
  two-simplex; carrier naturality alone leaves the full simplex, and \(9/5\)
  requires \(w_2=1/5\).
- The normalized Hilbert trace predicts \(17/9\), not \(9/5\).
- The \(4/5\) fractional law is enclosed on a continuous finite frequency
  window.
- Fractional wave evolution needs detector bandwidth for pointwise records.
- Three independent physical calibrations are necessary.
- Apparatus memory can be learned by a calibration-only nonparametric
  estimator in the synthetic validation.
- The local/fractional coupling obeys
  \(dx/d\log b=(6/5)x(1-x)\).
- Exact rational formulas reproducing \(\phi\) and \(\omega\) are rejected as
  target-driven arithmetic encodings.
- Projectivization removes \(\alpha_{\rm geo}\) from shape but neither
  completes the strict kernel nor transports amplitude-dependent roles.
- The damping factor is exact only conditional on the unresolved sector
  state.
- A signed receiver does not supply a strict signed preparation.
- No external physical dataset is admitted.

## Included files

- `FIN_Programs_125_137_Trace_Localizer_Physics_Monograph.pdf`
- `FIN_Programs_125_137_Trace_Localizer_Physics_Monograph.md`
- `FIN_Programs_125_137_Trace_Localizer_Physics_Monograph.tex`
- `FIN_Programs_125_137_Trace_Localizer_Physics_Results.json`
- `fin_programs_125_137_trace_localizer_physics.py`
- `test_fin_programs_125_137_trace_localizer_physics.py`
- `fin_programs_125_137_to_latex.py`
- `FIN_Programs_125_137_Trace_Localizer_Physics_Figures/`
- `FIN_Programs_125_137_Release_10_13_Manifest.json`

## Recommended Programs 138–150

The next research round should:

1. test modular/KMS selection of a sector state;
2. classify maximum-entropy reference measures;
3. test Morita stability of the exponent readout;
4. replace the guarded FFT by validated interval arithmetic;
5. derive a Diophantine oscillatory-tail remainder rate;
6. search for weighted fractional-wave dispersive estimates;
7. renormalize finite detector resolution;
8. prove joint system–instrument identifiability;
9. construct calibration-invariant observable ratios;
10. formulate a preparation resource theory for the sign torsor;
11. test a coupled state–damping variational principle;
12. build a typed full completion-map diagram;
13. freeze a pre-data physical protocol with rejection thresholds.

## Keywords

FIN; trace selection; cyclic groups; fibre localizer; stable process;
fractional Laplacian; wave propagation; ultraviolet cutoff; dimensional
calibration; apparatus tomography; information dynamics; kernel completion;
operational physics; mathematical physics.

## Suggested citation

Żuchowski, K. (2026). *FIN Programs 125–137: Trace Selection, Natural
Localizers, Fractional Physics, and Operational Sources* (FIN Research
Monograph, Release 10.13; Version 1.0.0) [Preprint]. Zenodo.
