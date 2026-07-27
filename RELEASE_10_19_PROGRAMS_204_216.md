# Release 10.19 — Categorical Measures, Perfect References, and External Falsification Gates

## Zenodo release description

**Title:** FIN Programs 204–216: Categorical Measures, Perfect References,
and External Falsification Gates

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-07-27  
**Language:** English  
**Publisher:** Zenodo  
**Access:** Open access  
**License:** Creative Commons Attribution 4.0 International (CC BY 4.0)

## Abstract

Release 10.19 reports thirteen executed research programs, numbered 204
through 216. The studies investigate categorical state selection,
target-independent exponent generation, formal proof infrastructure,
catalytic symmetry references, hidden dependence, anytime-valid open-set
testing, finite-shot process tomography, environmental moment bounds, phase
naturality, scale-free operator falsification, independent experimental
intake, and conditional prediction.

The release first proves a conditional central-state uniqueness theorem.
Inner-unitary naturality fixes normalized trace within every simple summand
of a finite semisimple algebra but leaves a central probability measure. If
one adds invariance under Morita autoequivalences that permute all simple
summands, the normalized central weights are uniquely uniform. For the FIN
dimension vector `(1,2,2,2,2)`, the corresponding exact expectation is
`9/5`. The normalized represented-space trace yields `17/9` and violates the
added permutation principle. This is not promoted to a strict state law:
the repository does not presently source Morita-permutation invariance as a
physical preparation rule.

A complementary no-go theorem shows that no automorphism-natural state
assignment is contravariantly natural under all unital `*`-homomorphisms. The
map

```text
f : C^2 -> C^3
f(a,b) = (a,b,b)
```

pulls the uniform state on `C^3` back to `(1/3,2/3)`, in conflict with the
uniform state `(1/2,1/2)` forced by automorphisms of `C^2`.

The strict exponent problem is classified within continuous one-parameter
cocycles. Every continuous solution of

```text
b(t+s) = b(t) + b(s)
```

has the form `b(t)=kappa*t`. Passing from legacy exponent `1` to strict
exponent `9/5` fixes only `kappa*t*=4/5`, leaving infinitely many rate-time
pairs. Multiplicative and stable-fixed-point flows have the same
nonselection defect. No target-independent strict exponent source or
legacy-to-strict bridge completion is obtained.

The strongest operational theorem concerns catalytic symmetry breaking. An
invariant catalyst cannot repair a violated one-copy reflection-conversion
trace-norm inequality. In contrast, the perfect reference

```text
Pi_+ = (I + sigma_x)/2
```

supports an explicit reflection-covariant measure-and-prepare channel that
maps the declared source to any prescribed target and returns `Pi_+`
exactly. Trace preservation, covariance, and mapping residuals are zero. The
construction is a sufficient operational selector, but the perfect
reference contains a fully asymmetric orientation bit. It supplies the
missing object and does not derive it from the strict kernel or discharge
`QW-2191`.

Further theorem-level contributions include:

- a Berbee-coupling plus DKW confidence theorem for hidden geometric mixing
  under a calibrated lower refresh bound;
- a fresh-block conformal e-process with exact time-uniform null control;
- a simultaneous finite-shot process-tomography region constructed from six
  binomial counts;
- a sharp environmental trigonometric moment theorem:
  `7/25 <= c2 <= 1` when `c1=4/5`, with explicit extremizing phase measures;
- an extension of the phase-source no-go from continuous circle
  homomorphisms to Borel homomorphisms and trivial-action continuous
  one-cocycles.

A scale-free spectral fingerprint is frozen before its internal challenge.
It divides the positive eigenvalues of the strict twelve-node generator by
the largest eigenvalue and uses a PSD/one-zero-mode structural gate. The
fingerprint is invariant across twelve decades of operator scaling and under
node permutation. It accepts a one-percent edge perturbation and rejects a
ten-percent perturbation, a nearest-neighbour cycle, and the raw legacy
generator. This is an internal operator falsification protocol, not an
experimental physical result.

The proof-infrastructure audit is explicit. A dependency-free Lean 4 core
probe compiles locally and 200 exact rational Dirichlet witnesses pass. The
general finite weighted-graph library does not compile because Mathlib is
absent from the clean local environment. Arb and python-flint are absent and
the Docker server is inaccessible. No formal sub-three-percent certificate
is claimed.

No local dataset passes the frozen independent eleven-field external intake
gate. The earlier synthetic method bundle remains synthetic and fails
independent provenance. Consequently the external prediction protocol is
frozen but correctly not executed.

## Main scientific results

1. **Morita-compatible central measure — Proven conditionally**
   - constraint rank: `4`;
   - normalized central weights: five copies of `1/5`;
   - exact selected expectation: `9/5`;
   - represented-space trace expectation: `17/9`;
   - all-homomorphism naturality defect: `1/6`;
   - strict state source: not obtained.

2. **Exponent-cocycle classification — Proven in class**
   - required increment: `4/5`;
   - displayed rate-time pairs: `601`;
   - maximum cocycle residual: `4.44e-16`;
   - target-independent `eta` source: not obtained.

3. **Arb environment — Contract complete; certificate open**
   - python-flint: absent;
   - Arb executable: absent;
   - Docker client: present;
   - Docker server: permission denied;
   - inherited worst enclosure: `0.6129831478`;
   - formal target below `0.03`: not certified.

4. **Lean proof layers — Partially machine checked**
   - local Lean version: `4.28.0`;
   - dependency-free core probe: compiled;
   - exact rational witness cases: `200`;
   - general Mathlib graph library: uncompiled because Mathlib is absent.

5. **Catalytic symmetry reference — Proven**
   - invariant-catalyst no-help theorem;
   - perfect-reference trace-preservation residual: `0`;
   - covariance residual: `0`;
   - exact mapping residual: `0`;
   - best sampled imperfect-catalyst necessary margin: `-0.808037`;
   - strict selector source: not obtained.

6. **Hidden-mixing inference — Proven under calibration**
   - record length: `10000`;
   - calibrated lower refresh probability: `0.05`;
   - selected thinning lag: `166`;
   - retained sample count: `61`;
   - nominal iid coverage: `0.845833`;
   - coupling-valid coverage: `1.0`;
   - rigorous mean interval width: `6.58742`.

7. **Sequential conformal e-process — Proven null law**
   - calibration records per step: `15`;
   - horizon: `20`;
   - Ville threshold: `20`;
   - iid Gaussian rejection: `0.02`;
   - sorted, AR(1), repeated-block, and drift rejection: `1.0`;
   - mean alternative crossing time: `2`.

8. **Finite-shot process region — Proven construction**
   - shots per phase and context: `1000`;
   - repetitions: `600`;
   - nonempty physical regions: `600`;
   - simultaneous empirical coverage: `1.0`;
   - memory-detection power: `1.0`.

9. **Environmental moment interval — Proven and sharp**
   - input moment: `c1=4/5`;
   - exact interval: `7/25 <= c2 <= 1`;
   - lower witness: equal atoms at `+/- arccos(4/5)`;
   - upper witness: mass `9/10` at `0` and `1/10` at `pi`.

10. **Phase naturality hierarchy — Proven no-go in class**
    - continuous homomorphisms: fail;
    - Borel homomorphisms: fail by automatic continuity;
    - trivial-action continuous one-cocycles: fail;
    - target-producing unrestricted interpolation: target-coded.

11. **Scale-free fingerprint — Frozen and challenged**
    - strict row sum: `1.6603072787660986`;
    - acceptance threshold: `0.02`;
    - all six exact scales and node permutation: accepted;
    - one-percent edge noise distance: `0.00233262`, accepted;
    - ten-percent edge noise distance: `0.0420993`, rejected;
    - nearest-neighbour cycle distance: `0.423325`, rejected;
    - raw legacy generator: rejected by PSD structural gate.

12. **Independent external intake — Gate failed**
    - admitted external bundles: `0`;
    - synthetic method bundle: `10/11`, not external;
    - web and Firecrawl: not used.

13. **Conditional external prediction — Correctly locked**
    - prerequisite external gate: false;
    - external prediction executed: no;
    - external physical validation claimed: no.

## Scientific interpretation

The same positive self-adjoint generator still supports wave, heat, and
Green dynamics:

```text
U_t = exp(-i t A)
P_t = exp(-t A)
G_z = (A + z I)^(-1).
```

Programs 204–216 show that this common engine does not uniquely generate the
objects required for an experiment. The centre of the observable algebra
requires a state-selection rule. Reflection-sensitive operations require an
asymmetry reference. Physical evolution requires a dimensional clock and
scale. Memory claims require interventions and calibrated acquisition.
Measurement requires an instrument and detector model. Validation requires
an independent ordered record.

The strongest interpretation surviving falsification is:

> FIN is a dimensionless spectral-information generator with an unresolved
> operational-section problem. Exact mathematics classifies parts of the
> required state, symmetry reference, environment, inference, and
> falsification structure, but the complete state–clock–instrument–record
> section is not forced by the strict operator.

## Scope guardrails

This release does **not** claim:

- a strict central state or temperature source;
- a target-independent source of `eta=9/5`;
- a completed legacy-to-strict bridge;
- transfer of legacy physical roles;
- a complete Arb certificate;
- complete Lean verification of the general graph library;
- a non-premise strict selector or discharge of `QW-2191`;
- a canonical physical unit;
- a role-bearing `L_total`;
- Standard Model or general-relativistic generation;
- external empirical validation;
- Theory-of-Everything closure.

## Recommended next studies

Programs 217–229 are specified in the monograph. The highest priorities are:

1. construct an explicit instrument-to-spectrum reconstruction theorem;
2. acquire one independent eleven-field external bundle;
3. quantify the cost and degradation of imperfect symmetry references;
4. extend the sharp environmental moment hierarchy;
5. complete the pinned Lean/Mathlib and Arb proof environments;
6. determine whether the Morita-permutation state principle has operational
   provenance;
7. execute the frozen external prediction only after the intake gate passes.

## Included files

- `FIN_Programs_204_216_Categorical_Catalytic_External_Monograph.pdf`
  — 48-page English preprint.
- `FIN_Programs_204_216_Categorical_Catalytic_External_Monograph.md`
  — editable source.
- `FIN_Programs_204_216_Categorical_Catalytic_External_Monograph.tex`
  — archival LaTeX.
- `FIN_Programs_204_216_Categorical_Catalytic_External_Results.json`
  — machine-readable results and claim boundaries.
- `fin_programs_204_216_categorical_catalytic_external.py`
  — executable research suite.
- `test_fin_programs_204_216_categorical_catalytic_external.py`
  — 77 regression and scientific-boundary tests.
- `FIN_Programs_204_216_Arb_Environment_Contract.json`
  — frozen certification environment.
- `FIN_Programs_204_216_Lean_Build_Contract.json`
  — formal build contract and local probe.
- `FIN_Programs_204_216_Lean_Core_Probe.lean`
  — locally machine-compiled dependency-free probe.
- `FIN_Programs_204_216_Dirichlet_Library.lean`
  — general Mathlib source, not yet compiled.
- `FIN_Programs_204_216_Phase_NoGo_Certificate.json`
  — naturality hierarchy certificate.
- `FIN_Programs_204_216_Scale_Free_Preregistration.json`
  — frozen projective spectral protocol.
- `FIN_Programs_204_216_External_Bundle_Request.json`
  — independent experimental intake contract.
- `FIN_Programs_204_216_External_Prediction_Preregistration.json`
  — locked held-out prediction protocol.
- `FIN_Programs_204_216_Categorical_Catalytic_External_Figures/`
  — thirteen figures.
- `FIN_Programs_204_216_Release_10_19_Manifest.json`
  — archival hashes.

## Reproduction

```bash
python3 fin_programs_204_216_categorical_catalytic_external.py
python3 -m unittest -v \
  test_fin_programs_204_216_categorical_catalytic_external.py
```

Expected result: 13 programs executed and 77 tests passing. The local
environment is expected to compile the dependency-free Lean core probe but
not the Mathlib-dependent general library unless Mathlib has been installed
and pinned.

## Keywords

spectral theory; operator algebra; Morita equivalence; categorical states;
Dirichlet forms; graph Laplacians; catalytic asymmetry; reference frames;
conformal e-processes; process tomography; trigonometric moment problems;
scale-free observables; operational physics; falsification.
