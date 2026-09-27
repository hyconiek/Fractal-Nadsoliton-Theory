# FINGERPRINT-OPTIMAL-PROBE-297
## A single lowest-harmonic scalar readout is sufficient for the declared D12 effective class, and a pre-hydrodynamic one-time protocol separates the report-293 FIN chain from the same-rho/same-exit comparator with controlled finite-sample error

Date: 2026-09-27

Status:
- exact sufficiency/tomography theorem inside the declared 12-state real symmetric D12-circulant effective class;
- exact finite-state likelihood calculations for the report-293 q_d inputs, N=3,...,8;
- explicit finite-sample Chernoff error bound under a declared readout-noise model;
- N=3,...,7 used for design, N=8 held out until the probe and time were frozen;
- not a laboratory validation and not a derivation of the q_d inputs from first principles.

## 1. Question

Report 293 showed that the resolved 12-state chain contains information beyond the coarse Z3 clock rho.  The remaining task is to find the lowest-complexity preparation/readout protocol that retains as much discrimination as possible between

1. the actual FIN effective shell rates q_1,...,q_6, and
2. the report-293 positive comparator with the same rho and the same total exit rate.

The comparator is reconstructed from the actual q by

    k_Z3 = q1+q2+q4+q5,
    rho  = 3 k_Z3,
    E    = 2(q1+q2+q3+q4+q5)+q6,

and

    q1^c=q2^c=q4^c=q5^c = k_Z3/4,
    q3^c=q6^c = (E-2 k_Z3)/3.

Equivalently, because 8(k_Z3/4)=2 k_Z3,

    E = 8 a + 3 b,
    a=k_Z3/4,
    b=(E-8a)/3.

This preserves exactly:
- rho;
- the coarse Z3 generator;
- the total 12-state exit rate.

It changes only the hidden shell decomposition.

## 2. Necessary no-go before optimization

The set of all positive q' with the same rho and total exit rate is a four-dimensional affine family locally around an interior q: six shell rates minus two independent scalar constraints.

Therefore alternatives q' satisfying the two calibration constraints can approach the FIN q arbitrarily closely.

By continuity of finite-time path/output distributions and of Chernoff information,

    inf_{q' != q, same rho, same exit} C(q,q') = 0.

Hence there is NO nonzero universal minimax error exponent against every distinct same-rho/same-exit countermodel unless a minimum alternative separation is imposed.

P1-297 is therefore optimized against the explicit report-293 equalized comparator, while the measurement theorem below is valid for the whole symmetric-circulant class.

## 3. Declared experimental noise model

Each independent shot is:

1. prepare the localized effective state J=0;
2. evolve for time t;
3. read one of the twelve localized states;
4. with probability 1-eta the label is read correctly;
5. with probability eta it is reported uniformly as one of the other eleven labels.

Thus

    P(Y_report=j | J_true=i)
      = 1-eta,        j=i,
      = eta/11,       j!=i.

The design value is

    eta = 0.05.

This 5% value is an illustrative declared noise level, not a measured apparatus specification.  Robustness is also reported for eta=0, 0.10 and 0.20.

The shots are assumed independent because each shot is separately prepared.  A single correlated long trajectory would require an effective-sample-size correction and is not covered by the quoted M bounds.

## 4. Exact one-mode sufficiency theorem

Let

    Y_1(J) = cos(2 pi J/12).

For J=0,...,11 this observable has exactly seven distinct values, corresponding to the reflection orbits

    {0}, {6}, {1,11}, {2,10}, {3,9}, {4,8}, {5,7}.

For a localized J=0 preparation and every real symmetric circulant generator,

    p_j(t)=p_{-j}(t).

The declared symmetric readout channel preserves the same reflection symmetry.

Therefore, for any two candidate generators in this class, the likelihood ratio is constant inside every pair {j,-j}.  Collapsing the twelve labels to the seven values of Y_1 loses no likelihood information.

Hence, for every t and eta in the declared model,

    C_full-label(t) = C_Y1(t),

and similarly the KL and all likelihood-ratio tests can be computed from Y_1 alone.

This is stronger than merely saying that k=1 has a large signal: it is a sufficient statistic for the complete one-time state readout in this symmetry class.

The k=5 cosine gives another one-to-one encoding of the same reflection orbits and is statistically equivalent.  k=1 is selected as the lowest harmonic.

### Consequence

No two-mode post-processing can improve the one-time discrimination information in this declared symmetric experiment, because Y_1 already attains the full-label information and data processing forbids an increase above it.

## 5. One-time shell tomography theorem

The same scalar observable does more than distinguish the report-293 comparator.

The seven-bin Y_1 histogram reconstructs the full symmetric transition row p_j(t), because each two-element orbit has equal probabilities.

For the symmetric circulant generator,

    lambda_k
      = 2 sum_{d=1}^5 q_d[cos(2 pi k d/12)-1]
        + q6[(-1)^k-1],

for k=1,...,6.

Starting from J=0,

    p_hat_k(t) = exp(lambda_k t).

Under the declared readout error,

    p_tilde_j = a p_j + eta/11,
    a = 1 - 12 eta/11.

For k!=0,

    p_tilde_hat_k = a exp(lambda_k t).

If eta is independently calibrated and a!=0, then

    lambda_k = t^(-1) log(p_tilde_hat_k/a).

Finally, the six lambda_k determine the six q_d because the 6x6 linear shell-to-spectrum matrix has

    rank = 6,

and exact determinant

    det M = -3456 sqrt(3) != 0.

Therefore, in the ideal infinite-sample version of the declared model:

    boxed:
    one localized preparation + one nonzero time + the seven-bin Y_1 histogram
    reconstructs the entire six-shell symmetric generator.

Finite data should use likelihood fitting rather than literal logarithmic inversion when empirical Fourier coefficients are noisy.

## 6. Time optimization without test reuse

To avoid using the same model instances for both tuning and evaluation:

- design set: N=3,4,5,6,7;
- held-out check: N=8.

The dimensionless time is

    tau = rho t.

For each candidate single cosine observable k=1,...,6, optimize the worst design-set Chernoff information

    C_k(tau)
      = -log min_{0<=s<=1}
          sum_y P_FIN(y;tau)^s P_cmp(y;tau)^(1-s).

The selected protocol is

    boxed:
    k=1,
    tau_* = 0.5427059873,
    t_* = 0.5427059873 / rho.

This is pre-hydrodynamic: the readout occurs at about 0.54 of the coarse Z3 relaxation time 1/rho.

The one-mode scan at eta=5% gives approximately:

| cosine mode | optimized worst-design Chernoff C |
|---|---:|
| k=1 | 0.00663362 |
| k=2 | 0.00210454 |
| k=3 | 0.00148485 |
| k=4 | 0 within numerical precision |
| k=5 | 0.00663362 |
| k=6 | 0.00008412 |

The zero information in k=4 is an internal negative control: k=4 is exactly the calibrated Z3 mode, so the two models were constructed to agree there.

## 7. Finite-sample error control

For equal prior probabilities and M independent shots, the standard Chernoff bound gives

    P_e^*(M) <= (1/2) exp(-M C).

Thus a sufficient count for the bound to fall below 5% is

    M >= ln(10)/C,

and below 1% is

    M >= ln(50)/C.

At the frozen tau_* and eta=5%:

| N | Chernoff C | M for bound <=5% | M for bound <=1% |
|---:|---:|---:|---:|
| 3 | 0.00675111 | 342 | 580 |
| 4 | 0.00663362 | 348 | 590 |
| 5 | 0.00683944 | 337 | 572 |
| 6 | 0.00727880 | 317 | 538 |
| 7 | 0.00792246 | 291 | 494 |

So the preregistered worst-design requirement is

    boxed: M=348 independent shots for a Chernoff upper bound <=5%.

## 8. Held-out N=8 result

After freezing k=1 and tau_*, N=8 was evaluated.

Result:

    C_N8 = 0.00871303018.

Therefore

    M >= 265

is sufficient for the same 5% Chernoff upper bound, and

    M >= 449

for 1%.

The seven Y_1 probabilities at N=8, eta=5%, tau=tau_* are:

| Y_1 | FIN | same-rho/same-exit comparator |
|---:|---:|---:|
| -1 | 0.0317665 | 0.0574896 |
| -sqrt(3)/2 | 0.0728324 | 0.0750878 |
| -1/2 | 0.1277624 | 0.0750878 |
| 0 | 0.1375985 | 0.1149792 |
| 1/2 | 0.0557753 | 0.0750878 |
| sqrt(3)/2 | 0.0439811 | 0.0750878 |
| 1 | 0.5302837 | 0.5271800 |

The discrimination does not come mainly from the probability of remaining at J=0.  It comes from the resolved pattern across several reflection orbits, which is exactly the hidden-shell information that the coarse Z3 observable erases.

## 9. Noise robustness with the probe frozen

Keeping the SAME k=1 and tau_* rather than retuning after seeing the noise level:

| eta | worst design C | held-out N=8 C | worst-design M for <=5% | N=8 M for <=5% |
|---:|---:|---:|---:|---:|
| 0% | 0.00803116 | 0.01042629 | 287 | 221 |
| 5% | 0.00663362 | 0.00871303 | 348 | 265 |
| 10% | 0.00548991 | 0.00728036 | 420 | 317 |
| 20% | 0.00374126 | 0.00503769 | 616 | 458 |

The fingerprint survives substantial symmetric readout noise; the cost is increased sample count.

## 10. Timing robustness

At eta=5%, keeping the probe fixed and perturbing tau by +/-20% gives worst-design Chernoff information

    tau=0.8 tau_*: C=0.00645497,
    tau=0.9 tau_*: C=0.00659214,
    tau=1.0 tau_*: C=0.00663362,
    tau=1.1 tau_*: C=0.00659779,
    tau=1.2 tau_*: C=0.00650026.

So the optimum is broad rather than knife-edge.  A moderate timing-calibration error does not destroy the test.

## 11. Operational protocol

A clean future test of the declared effective lane is therefore:

1. Use a separate calibration set to estimate rho from the coarse Z3 mode.
2. Pre-register the FIN q_d prediction and the comparator; do not estimate q_d from the discrimination shots.
3. Prepare localized state J=0 independently for every shot.
4. Evolve for

       t_* = 0.542706/rho.

5. Read J once.
6. Store only

       Y=cos(2 pi J/12),

   i.e. one of seven values.
7. Use the seven-bin likelihood-ratio test.
8. Calibrate eta independently or include it as a nuisance parameter in a preregistered likelihood model.

No second Fourier observable is required for one-time discrimination inside the declared symmetric-circulant class.

## 12. What 297 does and does not establish

### Established inside the declared effective model

- one scalar k=1 cosine readout is sufficient for the entire one-time symmetric state distribution;
- the same histogram is, in principle, enough to reconstruct all q_d;
- the report-293 FIN chain and equalized same-rho/same-exit comparator have a strictly positive finite-sample discrimination exponent;
- a concrete pre-hydrodynamic time and finite independent-shot budget are available;
- the result survives 10-20% symmetric readout corruption with larger sample counts;
- k=4 correctly carries zero comparator information because it is exactly the matched coarse mode.

### Not established

- the q_d are not laboratory-measured constants;
- the MZ reconstructed 12-state chain is still a conditional effective lane;
- no SI time scale is derived because rho is still an internal calibrated clock;
- no apparatus-specific eta is supplied;
- the quoted M counts assume independent preparations and exact candidate models;
- no finite sample can uniformly distinguish FIN from every arbitrarily close same-rho/same-exit alternative;
- one-time shell tomography does not test whether the 12-state dynamics is truly Markovian at multiple times.

## 13. New research implication

Report 297 changes the next experimental/mathematical question.

A second Fourier mode is not the next thing to add: k=1 already retains all one-time state information allowed by the declared symmetry.

The next genuinely new information must therefore come from MULTI-TIME structure.

A natural next task is:

    298 — MULTITIME-SEMIGROUP-FINGERPRINT

Test whether independently reconstructed Fourier factors satisfy

    p_hat_k(t2) = [p_hat_k(t1)]^(t2/t1)

or equivalently whether

    log p_hat_k(t)/t

is time-independent after correcting known readout noise.

Failure would detect residual Mori-Zwanzig memory/non-Markov structure that no one-time q_d fit can reveal.  This directly connects the fingerprint lane back to the memory cautions in reports 132-137, 259 and 296.

## Verdict

P1-297 is positive for the declared report-293 alternative and stronger than originally requested:

    boxed:
    the lowest cosine mode, used as a seven-outcome scalar observable,
    is a lossless one-time readout of the whole reflection-symmetric 12-state chain.

At eta=5%, a preregistered time

    tau_* = rho t_* = 0.5427059873

has worst-design Chernoff information 0.0066336191 and held-out N=8 information 0.0087130302.

The corresponding sufficient Chernoff-bound shot counts for <=5% equal-prior model-selection error are 348 in the worst design case and 265 in the held-out N=8 case.
