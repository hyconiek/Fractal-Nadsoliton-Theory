# 305 — PREPARATION-CONTRACT-PREDICTION-SET
## Set-valued microscopic FIN predictions remain separated from the same-rho comparator under an explicit preparation contract

Date: 2026-09-27

Status:
- exact finite-state `N=7` leave-one-out Gibbs propagation for the declared microscopic law;
- exact reuse of the D12-equivariant basin definition from reports 136/140 and tasks 300–304;
- convex-hull optimization over **all 2616 microscopic states in basin `J=0`** at the frozen late time;
- explicit controller family and preparation-resource accounting;
- finite-grid stress tests for clock/readout uncertainty;
- no asymptotic theorem and no laboratory apparatus calibration is claimed.

## 1. Why task 305 was necessary

Task 304 showed that the cross-N clock/shape law can predict the held-out `N=7`
seven-bin histogram to about 1.31% TV, while changing the microscopic
preparation inside the same macro-label `J=0` can move the prediction by more.
Therefore a statement of the form

    prepare J=0 -> evolve -> read Y

is not operationally complete until the microscopic preparation contract is
specified or its entire allowed prediction set is propagated.

This task does both.

Frozen microscopic point:

    N = 7
    g = 5.145228719489142

Frozen late protocol time from task 304's N<=6 clock law:

    rho_pred = 0.012502679779153265
    tau_design = 0.5427059873
    t_pred = 43.40717325296116

Nominal symmetric label-readout error:

    eta = 0.05.

The exact replay again gives:

    state count = 31824
    J=0 microscopic states = 2616
    total localized equilibrium mass = 0.9992530200759365
    equilibrium mass of J=0 = 0.08327108500632813
    stationarity residual = 7.70e-17.

## 2. Three canonical preparations are genuinely different microscopic contracts

All three have macro-label `J=0` at release:

1. **equilibrium-J0**
   `mu_eq(n)=pi(n)/pi(J0)` for `n in J0`;
2. **deep seed**
   all seven copies on label 0;
3. **flat-J0**
   uniform over the 2616 microscopic states classified into basin 0.

At the frozen late time, after exact microscopic FIN propagation and 5% readout
noise, their TV distances to the same-rho/same-total-exit comparator are:

    equilibrium-J0 : 0.0754337013
    deep seed      : 0.0826221841
    flat-J0        : 0.2966040747.

Thus the discriminator is not secretly relying on the deep seed; the broad
flat preparation actually moves farther from the declared comparator.

## 3. Operational pinning family

Introduce an explicit preparation controller that couples to the microscopic
occupation of label 0 while the system is restricted to the exact basin `J=0`:

    mu_kappa(n)
      proportional
    pi(n) exp[kappa n_0/N],        n in J0.

This is an **introduced operational controller model**, not a FIN-derived law.
It has useful endpoints:

    kappa=0       -> equilibrium Gibbs conditioned on J0,
    kappa->infty  -> deep seed.

A direct one-dimensional optimization gives the preparation in this entire
thermal-pin path closest to the comparator at the frozen late time:

    kappa_* = 0.85550762...
    TV_min  = 0.0730752580.

So the complete pinning path stays separated from the comparator by about
7.31% TV at nominal settings.

## 4. The strongest kill-test: arbitrary microscopic preparation inside J=0

Do not choose a preparation family at all.

For every pure microscopic state `i in J0`, compute its exact conditional
seven-bin output at the late time. Because propagation is linear, the set of
all predictions obtainable from **any probability distribution supported in
J0** is exactly the convex hull of these 2616 output vectors.

The minimum TV distance from the comparator to this whole convex set is a
linear program.

Result:

    boxed:
    min_{mu: supp(mu) subset J0}
        TV(P_FIN^mu(Y,t_pred), P_cmp(Y,t_pred))
      = 0.0632843086

at eta=0.05.

Therefore the comparator is **not** inside the full macro-only `J=0`
prediction set at the frozen nominal point.

The LP optimum is sparse, as expected for a convex problem. The nearest FIN
prediction uses approximately:

    91.7082% deep seed
     8.2918% state (3,0,2,0,0,0,0,0,0,0,2,0).

Even that optimized adversarial microscopic preparation cannot reproduce the
comparator.

This is much stronger than comparing equilibrium/seed/flat examples.

## 5. Preparation sensitivity first unfolds, then contracts

A subtle point is that preparation uncertainty is not monotonically erased
from `t=0`.

At `t=0` all three canonical preparations have the same observed basin label,
so the seven-bin diameter is zero. Once free evolution starts, hidden
microscopic differences become visible.

For the canonical set {equilibrium-J0, deep seed, flat-J0}, the maximum
pairwise TV diameter is approximately:

    t=0      : 0
    t=4      : 0.567636
    t=8      : 0.540982
    t=16     : 0.478254
    t=32     : 0.373558
    t=43.407 : 0.313645
    t=64     : 0.229400
    t=96     : 0.142071
    t=128    : 0.097234
    t=160    : 0.068408.

Thus the correct picture is:

    controller release
      -> hidden preparation information unfolds into observed transition data
      -> subsequent mixing contracts the preparation set.

The late initial-slip/preparation layer is therefore a dynamical object, not
merely a scalar nuisance.

Within the more physical pair equilibrium-J0 vs deep-seed the late difference
is much smaller:

    TV(t_pred) = 0.0222752

for the seven-bin conditional readout.

The large diameter is dominated by the strongly nonthermal flat preparation.

## 6. Explicit controller/resource accounting

Use relative entropy to the full equilibrium distribution as a dimensionless
preparation-resource measure:

    I_prep(mu) = D_KL(mu || pi)/ln 2  bits.

This is always an information cost. Under the additional standard thermal
control interpretation, `kT D_KL` is the reversible nonequilibrium free-energy
lower bound. FIN does **not** yet supply an SI value of `kT`, so no joule claim
is made.

Values at N=7:

    equilibrium-J0 : 3.58604057 bits
    deep seed      : 4.14385175 bits
    flat-J0        : 9.44317240 bits.

For the closest thermal-pin preparation:

    kappa_* ~= 0.8555
    I_prep ~= 3.5944 bits.

Hard selection of basin `J=0` itself has

    P_pi(J0) = 0.0832710850,
    surprisal = 3.58604057 bits,
    rejection-sampling mean attempts = 12.00897.

The additional microscopic sharpening from conditional equilibrium to the deep
seed costs only about 0.558 bits in this relative-entropy accounting, whereas
forcing a flat distribution over the whole basin is much more expensive.

This is the first explicit preparation-side resource ledger in this lane.
It does not yet include the physical implementation cost of measuring basin
membership, storing controller state, switching fields, or timing the release.
Those resources remain part of the autonomous-controller problem identified
after report 295.

## 7. Clock/readout stress test

For the explicit thermal-pin family, stress the protocol with:

    clock: t in [0.98,1.02] t_pred
    symmetric readout eta in [0,0.10].

These are declared stress ranges, not measured apparatus uncertainties.

On the sampled grid, using the **same calibrated nuisance parameters** for both
model hypotheses, the minimum direct-microscopic FIN/comparator separation is:

    0.0686041 TV.

Even if the two model prediction sets are conservatively allowed to choose
separate nuisance-grid values, the sampled minimum is:

    0.0661610 TV.

For the still broader arbitrary-`J0` convex set, exact LPs at the two clock
endpoints and midpoint give, before readout noise:

    0.0664867   at 0.98 t_pred
    0.0669353   at 1.00 t_pred
    0.0673904   at 1.02 t_pred.

At eta=0.10, common symmetric readout multiplies TV differences by

    1 - 12 eta / 11 = 0.890909...,

so the corresponding separations remain about 5.92--6.00% TV.

This is a **finite-grid finite-N certificate**, not a continuum-in-time theorem.

## 8. Important correction concerning the 304 cross-N error

Task 304 reported a 1.3069% TV cross-N prediction error at N=7 for the
**equilibrium-J0 preparation**.

It is tempting to subtract that number from every set-separation margin.
That would be too strong: task 304 did not prove a uniform 1.3069% error over
all microscopic preparations in `J=0`.

Therefore:

- direct microscopic FIN prediction sets are certified here;
- the explicit pin family is directly propagated and tested here;
- **uniform portability of the cross-N observable law over the whole
  preparation class remains open**.

This becomes the next mathematical target rather than being silently assumed.

## 9. Verdict

Task 305 is POSITIVE for the exact finite-N operational lane.

The strongest statement supported by the calculations is:

    exact microscopic FIN
      + declared macro preparation J=0
      + exact free evolution
      + declared seven-bin readout

produces a set of late-time predictions that remains separated from the
specified same-rho/same-total-exit comparator, even when the microscopic
preparation inside J=0 is optimized adversarially.

So the discriminator does **not** depend on silently selecting one microscopic
state inside `J=0`.

At the same time the task shows that preparation must stay in the theory:
its effects can be tens of percent over experimentally relevant windows, and
an external preparation controller consumes information/free-energy resources.

## 10. Epistemic boundary

This report does NOT prove:

- laboratory realizability of the controller;
- a physical meaning for labels, N, g or the microscopic clock;
- a uniform large-N preparation theorem;
- a uniform cross-N reduction bound over all preparations;
- that the chosen same-rho comparator exhausts alternatives;
- a fundamental law selecting the preparation controller;
- spatial composition, QM, GR or a fundamental FIN ontology.

The result is a controlled finite-state operational prediction-set theorem for
the declared mathematical models.

## 11. Next task

### 306 — PREPARATION-UNIFORM-CROSS-N-BOUND

The new bottleneck is no longer whether microscopic preparation can be hidden.
It is whether the compact cross-N observable law of task 304 approximates the
**entire preparation-conditioned FIN prediction set**, rather than only the
equilibrium-J0 center.

Target:

1. construct the effective initial-slip/preparation map from microscopic J0
   distributions into the 12-state effective description;
2. bound its error uniformly on a declared preparation class;
3. combine that bound with the cross-N clock/shape law;
4. prove or falsify disjointness of the resulting predicted set from the
   comparator set without direct microscopic propagation at the held-out N;
5. only then move to a finite-window/multi-time process theorem.
