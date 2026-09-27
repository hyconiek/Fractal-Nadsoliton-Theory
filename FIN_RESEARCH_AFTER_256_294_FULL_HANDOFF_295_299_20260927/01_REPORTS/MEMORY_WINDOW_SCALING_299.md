# MEMORY-WINDOW-SCALING-299
## The finite-N FIN memory window separates rapidly from the metastable Z3 clock, with an exact Stieltjes tail criterion

Date: 2026-09-27

Status:
- exact scalar Stieltjes-tail inequalities for the reversible Mori-Zwanzig memory sector;
- secondary finite-N analysis of the accepted N=3..8 slow-memory artifacts;
- independent coupled-hidden-gap check for N=6..8;
- descriptive finite-range scaling only: **no N→∞ theorem** and no certified asymptotic exponent.

Source Git blob SHAs:

- `crt_slow_memory_scaling_N3_N8.json`: `adc611c5da2bad0bbe73e7f3b56331cf58a1976d`
- `largeN_memory_N7.json`: `a36dd3dfb6258bf82a481bf96694c69910c08af0`
- `largeN_memory_N8.json`: `03b15853dd4cb807cd6ce525428f628a0b11c28d`
- `rigorous_memory_tail_bounds_N6.json`: `3c5923466d1da889e5e9b1ff602f98adb75cde1b`
- `STIELTJES_MEMORY_STRUCTURE_231.md`: `15d3f314c7fe439dee76fa64950cccef3c20e987`

## 1. Question after report 298

Report 298 showed at N=6 that the exact projected microscopic correlation has a short multi-rate boundary layer and then converges extremely accurately to one slow exponential / one effective D12 generator.

The next question is whether the width of that memory layer remains relevant when the metastable clock becomes slower with N.

The correct dimensionless quantity is not the microscopic memory time by itself. It is

\[
\boxed{\epsilon_{\rm mem}(N)=\rho_N\,\tau_{\rm mem}(N)}
\]

where

\[
\tau_{\rm mem}=\frac{M_1}{M_0}
\]

and \(\rho_N=-\lambda^{\rm MZ}_{k=4}\) is the exact Z3 clock of the accepted 12-state effective generator.

If \(\epsilon_{\rm mem}\ll1\), memory occupies only a small fraction of the slow metastable time \(1/\rho\).

## 2. Exact Stieltjes tail theorem

Report 231 proved for one reversible symmetry sector

\[
K(t)=\int_0^\infty e^{-\gamma t}\,d\mu(\gamma),\qquad d\mu\ge0.
\]

Define

\[
M_0=\int_0^\infty K(t)dt,
\qquad
M_1=\int_0^\infty tK(t)dt.
\]

Then

\[
f(t)=\frac{K(t)}{M_0}
\]

is a probability density on \(t\ge0\), with mean

\[
\mathbb E_f[t]=\frac{M_1}{M_0}=\tau_{\rm mem}.
\]

Therefore Markov's inequality gives the exact integrated-tail bound

\[
\boxed{
\frac{1}{M_0}\int_T^\infty K(t)dt
\le
\frac{\tau_{\rm mem}}{T}.
}
\]

Put

\[
T=\frac{\theta}{\rho},
\]

i.e. wait a fraction \(\theta\) of one metastable Z3 relaxation time. Then

\[
\boxed{
\mathrm{Tail}(\theta/\rho)
\le
\frac{\epsilon_{\rm mem}}{\theta},
\qquad
\epsilon_{\rm mem}=\rho\frac{M_1}{M_0}.
}
\]

Thus \(\epsilon_{\rm mem}\to0\), if eventually proved, would be a direct sufficient criterion that the integrated memory boundary layer becomes negligible on every fixed positive fraction of the slow time.

This implication is exact. What remains unproved is the asymptotic premise \(\epsilon_{\rm mem}\to0\).

## 3. Finite-N result: memory time stays microscopic while the slow clock expands

Accepted values in the k=4 / Z3 memory sector are:

| N | rho | 1/rho | M1/M0 | epsilon_mem=rho M1/M0 |
|---:|---:|---:|---:|---:|
| 3 | 0.1320059567 | 7.5754 | 0.371875 | 0.0490898 |
| 4 | 0.0718543830 | 13.9170 | 0.359667 | 0.0258437 |
| 5 | 0.0402247567 | 24.8603 | 0.354449 | 0.0142576 |
| 6 | 0.0226189122 | 44.2108 | 0.309489 | 0.00700031 |
| 7 | 0.0126419111 | 79.1020 | 0.306945 | 0.00388037 |
| 8 | 0.00698992203 | 143.063 | 0.270824 | 0.00189304 |

Across N=3..8:

- the slow time \(1/\rho\) grows by a factor **18.89**;
- the mean memory time changes only from 0.372 to 0.271, a factor **1.37**;
- \(\epsilon_{\rm mem}\) falls by a factor **25.93**.

The same computation using the independently extracted exact microscopic slow eigenvalue rather than \(\rho_{\rm MZ}\) gives essentially identical values; at N=8 it is 0.00189297.

So over the entire tested range, the decrease is not caused by defining the clock through the MZ approximation.

## 4. Exact integrated-tail consequences

Using only the moment theorem, after waiting 10% of one slow time,

\[
T=0.1/\rho,
\]

the integrated unresolved memory tail is bounded by:

| N | Tail bound after 0.1/rho |
|---:|---:|
| 3 | 49.1% |
| 4 | 25.8% |
| 5 | 14.3% |
| 6 | 7.00% |
| 7 | 3.88% |
| 8 | **1.89%** |

At N=8, after 20% of a slow time, the moment-only bound is already below **0.947%**.

These are conservative inequalities. They use only M0 and M1 and make no lower-gap assumption.

## 5. Stronger coupled-gap certificate at N=6..8

If every hidden mode that actually couples to this memory sector obeys

\[
\gamma\ge\gamma_c>0,
\]

then positivity of the Stieltjes measure gives the stronger exact inequality

\[
\frac{1}{M_0}\int_T^\infty K(t)dt
\le e^{-\gamma_cT}.
\]

Equivalently, define

\[
\epsilon_{\rm gap}=\frac{\rho}{\gamma_c}.
\]

Then at \(T=\theta/\rho\),

\[
\boxed{
\mathrm{Tail}\le e^{-\theta/\epsilon_{\rm gap}}.
}
\]

The accepted finite-state values are:

| N | gamma_c | rho/gamma_c |
|---:|---:|---:|
| 6 | 0.677188 | 0.0334012 |
| 7 | 0.484279 | 0.0261046 |
| 8 | 0.499753 | **0.0139868** |

Thus at 10% of one slow time the exponential tail bounds are:

- N=6: **5.01%**;
- N=7: **2.17%**;
- N=8: **0.0785%**.

At N=8, 20% of one slow time gives

\[
\boxed{6.16\times10^{-7}}
\]

for the remaining normalized integrated memory tail.

Equivalently, using this finite-N coupled gap, the waiting fraction needed to force the tail below 1% is at most:

- N=6: 0.154 of one slow time;
- N=7: 0.120;
- N=8: **0.0644**.

For a 0.1% tail the N=8 bound requires only **0.0966** of one slow time.

Caveat: the N=6 gap comes from the stored proof-grade tail-bound calculation; the N=7 and N=8 coupled-gap values are accepted finite-state numerical eigensolve outputs, not interval-certified large-N theorems.

## 6. Initial slip and MZ error shrink independently

Two other quantities improve monotonically over the same N=3..8 range.

The one-pole Stieltjes residue is

\[
Z=\frac{1}{1+M_1}.
\]

The initial-slip deficit \(1-Z\) falls:

\[
0.0811,\ 0.0717,\ 0.0398,\ 0.0273,\ 0.0146,\ 0.00858.
\]

The relative error of the M0+M1 slow eigenvalue versus the exact microscopic slow eigenvalue falls:

\[
4.31\times10^{-3}
\to
3.58\times10^{-5}.
\]

So three different diagnostics improve together:

1. memory duration relative to the slow clock;
2. initial-slip amplitude;
3. error of the local MZ Markov closure.

This makes the scale-separation interpretation substantially harder to attribute to one arbitrary diagnostic.

## 7. Descriptive scaling, not an asymptotic law

A log-linear fit over only N=3..8 gives

\[
\epsilon_{\rm mem}(N)
\approx
0.3486\,e^{-0.6479N},
\]

with log-space \(R^2\approx0.9993\), corresponding to a factor about 0.523 per added copy.

The MZ relative error has a descriptive factor about 0.381 per added copy.

These fits are useful for planning the next computation only. Six finite-N points do **not** establish an exponential asymptotic law, an Eyring-Kramers exponent, or \(\epsilon_{\rm mem}\to0\).

## 8. Scientific verdict

Report 299 is positive in the finite-N sense:

\[
\boxed{
\text{for N=3..8, microscopic memory remains O(1) while the effective Z3 clock becomes rapidly slower.}
}
\]

More sharply,

\[
\boxed{
\epsilon_{\rm mem}=\rho M_1/M_0
\text{ decreases monotonically from }4.91\%\text{ to }0.189\%.
}
\]

The coupled-gap calculation at N=6..8 independently supports the same separation.

Therefore the emerging architecture is now:

\[
\text{exact finite-N Gibbs dynamics}
\to
\text{short reversible memory layer}
\to
\text{controlled localized 12-state generator}
\to
\text{exact Z3 quotient}.
\]

This is stronger than merely observing that the memory kernel decays: the relative width of the memory layer is demonstrably shrinking against the endogenous metastable clock over the tested sequence.

## 9. What is still missing

The decisive missing theorem is now very specific.

One must show either

\[
\inf_N\gamma_c(N)>0
\]

for the hidden modes that actually couple to the Z3 memory sector while \(\rho_N\to0\), or directly prove

\[
\boxed{\rho_N M_1(N)/M_0(N)\to0.}
\]

The global microscopic spectral gap is not the correct object because symmetry-decoupled metastable modes can approach zero without carrying memory weight.

The correct next target is therefore:

**300 — COUPLED-MEMORY-GAP-LARGEN:** extend the coupled-sector computation to N=9..12 if feasible and search for an N-uniform lower bound on the symmetry-allowed hidden memory spectrum.

A successful lower bound together with the already observed metastable collapse of \(\rho_N\) would turn the present finite-N pattern into an actual asymptotic Markov-separation theorem.

## 10. Scope boundary

Nothing in report 299:

- proves the N→∞ limit;
- derives a physical spatial site or dimension;
- supplies SI time or an experimental apparatus;
- solves the complete-system source problem from report 295;
- derives QW-2191, the legacy-to-strict bridge, Standard Model, gravity, or a ToE.

It is a controlled finite-N dynamical result inside the declared leave-one-out Gibbs / localized-basin lane.
