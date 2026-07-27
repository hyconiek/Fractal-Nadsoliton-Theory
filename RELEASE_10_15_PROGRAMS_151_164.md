# Release 10.15 — FIN Programs 151–164

## Axiomatic Operational Foundations, State Selection, and Falsifiable Measurement

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

Release 10.15 executes fourteen mathematical, computational, and operational
research programs and reorganizes the FIN framework axiomatically.

The minimal common mathematical object is identified as a conservative
spectral Dirichlet generator. Its functional calculus produces the unitary
and diffusive dynamics

\[
U_t=e^{-itA},
\qquad
P_t=e^{-tA},
\]

together with the Green family. The same generator does not, however,
uniquely determine a physical state, signed preparation, measurement
instrument, dimensional calibration, or legacy-to-strict completion map.
Explicit countermodels prove that these are independent obligations.

The release constructs one economical new state-selection axiom. For fibre
dimensions \((1,2,2,2,2)\), take

\[
H_F|_{\mathcal H_p}=(d_p-1)I
\]

and

\[
A_{\rm ME}:\qquad
\beta_F
=
\frac{\alpha_{\rm geo}}
{|\operatorname{Hom}(U(12),\{\pm1\})|}
=\frac{\alpha_{\rm geo}}4
=\log2.
\]

The Hilbert-degenerate Gibbs state then has exactly uniform central weights:

\[
w_p
=
\frac{d_p e^{-\beta_F(d_p-1)}}
{\sum_q d_q e^{-\beta_F(d_q-1)}}
=\frac15,
\qquad
\eta=\sum_p w_pd_p=\frac95.
\]

This conclusion is unique for the declared unit energy gap, but remains
explicitly axiom-augmented rather than strict-derived.

A six-axiom architecture, AFIS, is then constructed:

1. conservative spectral Dirichlet generator;
2. sector equilibrium/state;
3. preparation resource;
4. instrument, environment, apparatus, and record;
5. dimensional calibration;
6. typed legacy-to-strict completion.

All \(2^6=64\) subsets are evaluated mechanically. Removing every axiom has a
specific countermodel and destroys a named capability. The smallest calibrated
signed test requires the first five axioms; the FIN-specific completion remains
a separate sixth obligation.

## Principal results

- A formal \(N=32768\) interval FFT covers 100 continuous frequency cells.
  Every symbol lower bound is positive and every cell is compatible with
  \(C|q|^{4/5}\). The worst formal remainder improves to \(0.612983\), but a
  sub-\(3\%\) theorem is not yet obtained.
- Irrationality alone is proved insufficient for an all-scale polynomial
  discrepancy rate. The missing premise is isolated as a Diophantine
  condition \(\|q\theta\|\ge\kappa q^{-\nu}\).
- Invariant fibre-groupoid states form a two-dimensional simplex:
  \((x,y,z/3,z/3,z/3)\). Naturality does not select uniform sector weights.
- \(A_{\rm ME}\) selects \(w_p=1/5\) and \(\eta=9/5\) exactly and uniquely in
  the declared Gibbs model.
- The reflection monotone \(M(\rho_r)=|r|\) is complete for deterministic
  conversions on the two-branch orbit line.
- The stable-IQR asymptotic variance is derived and validated numerically.
  At \(n=4000\), the predicted slope SD is \(0.0307072\) and the Monte Carlo
  SD is \(0.0294257\).
- An unrestricted multiplicative apparatus response is proved to make the
  spreading exponent nonidentifiable.
- The frozen Program-150 threshold is adversarially falsified as a
  FIN-specific test: generic stable scaling and time-dependent apparatus gain
  can trigger it.
- The legacy period-eight character and the strict infinite-order phase are
  proved not to admit a generator-preserving equivariant character bridge.
- Exhaustion of 343 structural energy formulas makes \(d-1\) the unique
  shortest target-realizing formula in the declared grammar.
- A six-edge obstruction matrix licenses zero current transfers of nine
  audited legacy roles.
- The external-data schema is executable, but zero datasets are admitted.
- The AFIS independence theorem proves that no single theorem acting on the
  bare spectral generator can create all missing physical structures.

## Executed programs

1. **Program 151:** tighter validated fractional certificate.
2. **Program 152:** axiomatic Diophantine condition and irrationality-only
   obstruction.
3. **Program 153:** functorial fibre-groupoid probability classification.
4. **Program 154:** axiom-augmented modular equilibrium.
5. **Program 155:** completeness of the reflection resource monotone.
6. **Program 156:** detector deconvolution and finite-resolution bias.
7. **Program 157:** semiparametric system–instrument identifiability.
8. **Program 158:** finite-sample stable-IQR inference.
9. **Program 159:** blind adversarial protocol challenge.
10. **Program 160:** phase/frequency bridge obstruction theorem.
11. **Program 161:** exhaustive reference-energy grammar.
12. **Program 162:** conditional legacy-role obstruction matrix.
13. **Program 163:** external intake readiness.
14. **Program 164:** minimal axiomatic system and independence countermodels.

## Scientific interpretation

The release supports an axiomatic informational-spectral operational
interpretation of FIN:

\[
\text{spectral information}
\to
\text{selected state}
\to
\text{preparation and instrument}
\to
\text{calibrated record}.
\]

The spectral generator unifies mathematical propagation and diffusion.
Physical prediction begins only after the independent operational and
dimensional layers are supplied. The legacy-to-strict completion remains
FIN-specific and incomplete.

A mirror/reflection coupling represents the relevant \(\mathbb Z_2\)
symmetry, but a symmetric mirror pair does not choose its own sign. A signed
preparation, inversion-odd source, noncovariant instrument, or explicit
boundary condition remains necessary. QW-2191 is not discharged.

## Included files

- `FIN_Programs_151_164_Axiomatic_Operational_Foundations_Monograph.pdf`
- `FIN_Programs_151_164_Axiomatic_Operational_Foundations_Monograph.md`
- `FIN_Programs_151_164_Axiomatic_Operational_Foundations_Monograph.tex`
- `FIN_Programs_151_164_Axiomatic_Operational_Results.json`
- `FIN_Programs_151_164_Minimal_Axiomatic_System.json`
- `fin_programs_151_164_axiomatic_operational_foundations.py`
- `test_fin_programs_151_164_axiomatic_operational_foundations.py`
- `fin_programs_151_164_to_latex.py`
- `FIN_Programs_151_164_Axiomatic_Operational_Figures/`
- `FIN_Programs_151_164_Release_10_15_Manifest.json`

## Recommended Programs 165–177

The next round should:

1. certify the oscillatory tail with Arb and cancellation-aware bounds;
2. derive or falsify \(A_{\rm ME}\) from composition and coarse graining;
3. formalize AFIS consistency and independence in a proof assistant;
4. seek an effective all-scale irrationality measure;
5. classify monoidal fibre-groupoid valuations;
6. prove perturbative stability of the \(A_{\rm ME}\) equilibrium;
7. extend reflection resource theory to the full density space;
8. derive nonasymptotic detector statistics;
9. optimize calibration-control experiments;
10. preregister a composite adversarial model comparison;
11. generalize the categorical phase obstruction;
12. test exactly one explicit nonlinear completion candidate;
13. build an AFIS CP-instrument model of double-slit records.

## Guardrail

This release does not claim:

- strict derivation of \(A_{\rm ME}\);
- strict selector closure or QW-2191 discharge;
- internal generation of physical units;
- completion of the legacy-to-strict bridge;
- transfer of any legacy physical role;
- a role-bearing \(L_{\rm total}\);
- Standard Model or general-relativistic generation;
- Theory-of-Everything closure;
- external experimental validation.

## Keywords

FIN; axiomatic mathematical physics; Dirichlet forms; spectral theorem;
operator semigroups; modular equilibrium; KMS states; information theory;
resource theory; quantum instruments; dimensional calibration; fractional
dynamics; semiparametric identifiability; detector physics; kernel completion.

## Suggested citation

Żuchowski, K. (2026). *FIN Programs 151–164: Axiomatic Operational
Foundations, State Selection, and Falsifiable Measurement* (FIN Research
Monograph, Release 10.15; Version 1.0.0) [Preprint]. Zenodo.
