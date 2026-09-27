# MICROSCOPIC-PERIODIC-BRIDGE-296
## Exact finite cyclic bridge for the finite-N leave-one-out Gibbs process, before the metastable 12-state/Z3 reduction

Date: 2026-09-27

Status:
- exact finite-state bridge theorem applied directly to the accepted finite-N labelled Gibbs heat-bath process;
- exact uniformization construction and exact cylinder formula;
- explicit spectral cylinder-error bound;
- direct full-labelled numerical replay for N=2,3 and exact occupation-count quotient replay through N=7;
- no claim that the periodic history carrier is the physically complete simultaneous carrier.

This report implements the next task after COMPLETE-SYSTEM-CRITERION-295.
It deliberately does **not** insert the metastable 12-state chain or its Z3 quotient.

---

## 1. Microscopic process used

Use the exact finite-N labelled-copy Gibbs law already admitted in report 61:

    Pi_N(x_1,...,x_N)
      proportional to
      12^(-N)
      exp[(g/(2N)) sum_(a,b) A7[x_a,x_b]].

For one chosen label a, let m be the occupation counts of the other N-1 labels.
The exact leave-one-out conditional is

    q_j^(-a) = softmax_j[(g/N)(A7 m)_j].

Each label has rate 1 and, when selected, is redrawn from this exact conditional.
The continuous-time generator Q_N is reversible with respect to Pi_N.

For the controlled replay below use the already occurring coexistence/balance point

    g = 5.145228719489144.

Nothing in the analytic bridge theorem depends on this particular g.

---

## 2. Clean discrete-time skeleton

Let K_N be the one-event kernel:
- choose one of the N labels uniformly;
- redraw it from its exact leave-one-out Gibbs conditional.

Then

    Q_N = N (K_N - I).

Uniformize at rate 2N:

    boxed:
    P_N = I + Q_N/(2N) = (I + K_N)/2.

This is useful for two reasons.

1. It is an exact skeleton of the same continuous-time chain.
2. Since the exit rate of Q_N is at most N, all eigenvalues of Q_N lie in [-2N,0], so P_N has spectrum in [0,1].

Thus P_N is primitive, reversible and spectrally nonnegative.
No coarse Q3 law has been introduced.

---

## 3. Finite periodic microscopic bridge

For a period L define a probability measure on L complete microscopic states

    x_0,...,x_(L-1)

by

    boxed:
    mu_(N,L)(x_0,...,x_(L-1))
      = [1/Z_(N,L)] product_(t=0)^(L-1) P_N(x_t,x_(t+1)),

with

    x_L=x_0

and

    Z_(N,L)=tr(P_N^L).

The deterministic cyclic shift

    (x_0,x_1,...,x_(L-1))
      ->
    (x_1,...,x_(L-1),x_0)

is a bijection preserving mu_(N,L).

Therefore the finite bridge is an **exact finite reversible carrier** of the microscopic path process.

This is the microscopic analogue of report 290, not a construction made after metastable reduction.

---

## 4. Exact cylinder formula

Let the stationary two-sided microscopic process have cylinder law

    mu_inf(x_0,...,x_r)
      = Pi_N(x_0) product_(t<r) P_N(x_t,x_(t+1)).

For r<L, summing the periodic bridge over the unobserved part of the cycle gives

    boxed:
    mu_(N,L)(x_0,...,x_r)
      = [product_(t<r) P_N(x_t,x_(t+1))]
        (P_N^(L-r))_(x_r,x_0)
        / tr(P_N^L).

Hence

    boxed:
    mu_(N,L) / mu_inf
      =
    (P_N^(L-r))_(x_r,x_0)
    /
    [Pi_N(x_0) tr(P_N^L)].

This is exact for every finite N, every L>r and every allowed cylinder.

---

## 5. Explicit convergence bound

Because P_N is finite, reversible and nonnegative-spectrum, choose an orthonormal L2(Pi_N) eigenbasis with

    1=lambda_1 > lambda_2 >= ... >= lambda_M >= 0.

Write

    lambda_* = lambda_2 < 1.

The reversible spectral expansion gives

    P_N^m(y,x)/Pi_N(x)
      =
    1 + sum_(k>=2) lambda_k^m psi_k(y) psi_k(x).

By Cauchy-Schwarz,

    |P_N^m(y,x)/Pi_N(x)-1|
      <=
    lambda_*^m
    sqrt[(1/Pi_N(x)-1)(1/Pi_N(y)-1)].

Also

    tr(P_N^L)=1+sum_(k>=2) lambda_k^L,

so

    0 <= tr(P_N^L)-1 <= (M-1) lambda_*^L.

Combining these with the exact cylinder ratio yields the explicit pointwise relative-error bound

    boxed:
    |mu_(N,L)/mu_inf - 1|
      <=
    lambda_*^(L-r)
    sqrt[(1/Pi_N(x_0)-1)(1/Pi_N(x_r)-1)]
    + (M-1) lambda_*^L.

A cruder uniform bound follows by replacing the square root by

    1/Pi_min - 1.

Therefore every fixed microscopic cylinder converges to the stationary two-sided microscopic process as L->infinity.

This satisfies the requested bridge criterion without using Q12 or Q3.

---

## 6. Separation of bridge convergence from metastable mixing

Let Delta_N be the continuous-time spectral gap of -Q_N.
Since

    P_N = I + Q_N/(2N),

we have

    boxed:
    lambda_* = 1 - Delta_N/(2N).

Therefore

    lambda_*^(L-r)
      <=
    exp[-Delta_N (L-r)/(2N)].

If one uniformized step represents mean microscopic time 1/(2N), define

    T_close=(L-r)/(2N).

Then the closing error decays as

    exp[-Delta_N T_close].

This gives the requested clean separation:

- the periodic-bridge construction itself is exact for every finite N;
- the **period needed for local convergence** is controlled by the microscopic relaxation/mixing scale 1/Delta_N;
- metastability enters only through a small Delta_N and can make the required period very large.

Thus one must not say merely "L is much larger than the observed window".
In a metastable regime the correct requirement is roughly

    T_close >> 1/Delta_N,

up to the desired logarithmic accuracy and endpoint-weight factors.

---

## 7. Noninteracting control: uniformization is not causing the slowdown

At g=0 every label independently refreshes to the uniform 12-state law at rate 1.
The labelled product chain then has exact continuous-time spectral gap

    Delta_N(g=0)=1

for every N.

The lazy skeleton gap is therefore

    gap(P_N)=1/(2N),

but after conversion back to microscopic time the relaxation time is exactly 1.

Hence any large growth of 1/Delta_N at the interacting g used below is not an artefact of choosing the 2N-uniformized skeleton.

---

## 8. Full labelled microscopic replay

The supplied script constructs the **full labelled chain**, not only occupation counts.

At

    g = 5.145228719489144

it gives:

### N=2

    labelled states          = 144
    gap(P_N)                 = 0.06472285452741677
    Delta_N                  = 0.2588914181096671
    continuous tau=1/Delta   = 3.862623208222365
    max detailed-balance residual
                              = 6.51e-19

### N=3

    labelled states          = 1728
    gap(P_N)                 = 0.02115761613746803
    Delta_N                  = 0.1269456968248082
    continuous tau=1/Delta   = 7.877383991834346
    max detailed-balance residual
                              = 2.71e-19

These are floating finite-matrix replays, not interval eigenvalue certificates.

For N=2 a direct dense evaluation of the complete cylinder TV distance gives, for a 3-state observed block (r=2):

    L=8      TV = 1.7487053905e-1
    L=32     TV = 9.5546477139e-2
    L=64     TV = 1.5316433490e-2
    L=128    TV = 1.7164098838e-4
    L=256    TV = 2.4653346551e-8
    L=512    TV = 1.24e-15  (floating floor)

The finite bridge therefore approaches the true microscopic two-sided process exactly as predicted by the spectral theorem.

The TV sequence need not be monotone at short periods because the normalization tr(P^L) and multiple modes contribute; the theorem concerns the controlled large-L decay.

---

## 9. Exact count quotient used only as a larger-N spectral diagnostic

The occupation-count chain is an exact Markov quotient of the labelled Gibbs process, not the later metastable 12-state reduction.
It was used to continue the spectral diagnostic to N=7.

At the same g:

| N | count states | Delta_N | tau=1/Delta_N |
|---:|---:|---:|---:|
| 2 | 78 | 0.2588914181 | 3.8626 |
| 3 | 364 | 0.1269456968 | 7.8774 |
| 4 | 1365 | 0.06862964088 | 14.5710 |
| 5 | 4368 | 0.03817290434 | 26.1966 |
| 6 | 12376 | 0.02132466145 | 46.8941 |
| 7 | 31824 | 0.01184068929 | 84.4545 |

For N=2 and N=3 the slow gap of the full labelled chain agrees numerically with the count quotient to <2e-11.
Do not extrapolate that equality to all N without proof.

The N=2..7 values show a strong finite-N growth of the relaxation time near this interacting point.
They are consistent with an emerging metastable slowdown, but **do not prove an asymptotic exponential law or B4 exponent**.

This is exactly why the periodic length and metastable mixing time must be kept separate.

---

## 10. Direct connection to report 295

Condition on which microscopic label is updated.
The current full labelled microstate still has stochastic branching because the new symbol is sampled from q^(-a).

The stationary average fresh-output entropy is numerically

    N=2: 1.759730540691066 bits / conditioned update
    N=3: 1.302008749716988 bits / conditioned update.

Thus the present microstate alone has

    F_e = H(X'|X,e) > 0.

By report 295 it is not a complete reversible state carrier.

In contrast, on the full periodic path carrier the elementary shift is a deterministic permutation, so

    F_shift = B_shift = 0.

Therefore 296 gives an explicit microscopic realization of the 295 logic:

    stochastic current state
      -> hidden path information is required;

    complete cyclic history carrier
      -> reversible zero-defect shift.

But this also sharpens the no-go:

    the history carrier can always absorb the apparent fresh randomness.

So the bridge does **not** prove that those history slots are physical simultaneous subsystems.

---

## 11. Does information continuity survive before Q3 reduction?

Yes, in a precise mathematical sense.

The chain

    exact finite-N labelled Gibbs law
      -> exact lazy microscopic skeleton P_N
      -> stationary two-sided microscopic path process
      -> finite cyclic bridge
      -> deterministic reversible shift

is valid without:
- Q12;
- Z3;
- metastable basin labels;
- a supplied coarse reset kernel.

Therefore the reversible information-continuity construction is not an artefact of the metastable Q3 reduction.

However what survives is a **path-space reversible completion theorem**, not yet a physical locality theorem.

---

## 12. What 296 does NOT establish

296 does not show that:
- the L history positions are simultaneous physical sites;
- FIN itself selects a finite L;
- the cyclic history shift is the fundamental microscopic law;
- the finite-N Gibbs heat-bath randomness is fundamentally generated by those history slots;
- a spatial dimension has been derived;
- SWAP is compulsory on physical simultaneous units;
- the observed finite-N slowdown already proves a large-N metastable exponent;
- QW-2191, SI units, the legacy-to-strict bridge, L_total, SM, GR or a ToE are closed.

The construction is a mathematically exact reversible dilation of the accepted microscopic stochastic process.
Report 295 still forbids promoting that representation to an absolutely complete physical carrier without an independent FIN source of the carrier boundary.

---

## 13. Verdict

**MICROSCOPIC-PERIODIC-BRIDGE-296 passes.**

The important result is stronger than simply repeating report 290:

    boxed:
    the finite cyclic/natural-extension information-continuity picture already exists at the exact finite-N leave-one-out Gibbs level, before any metastable Q12/Z3 coarse graining.

The exact convergence scale is controlled by the microscopic spectral gap:

    boxed:
    bridge error ~ exp[-Delta_N T_close].

So metastability does not invalidate the bridge, but can make the finite period required for a faithful local approximation parametrically large.

At the same time 296 reinforces the central limitation of 295:

    reversible history completion is mathematically available,
    but does not identify history coordinates with physical simultaneous degrees of freedom.

---

## 14. Next research move

The handoff's next numbered task is

    FINGERPRINT-OPTIMAL-PROBE-297.

That task can now be sharpened:

1. use the actual q_d(N) / resolved effective chain;
2. construct same-rho and same-total-exit countermodels;
3. choose the smallest one- or two-mode preparation/readout protocol;
4. include declared measurement noise;
5. optimize discrimination in a pre-hydrodynamic time window;
6. keep training/tuning data separate from the final discrimination test.

In parallel, the foundational lane exposed by 295–296 remains:

    AMBIENT-SIMULTANEOUS-CARRIER-SOURCE:
    what FIN-internal object distinguishes a physical simultaneous carrier from a reversible history dilation?
