# STRENGTHENED-RECURSIVE-LOCAL-CHILD-TEST-131
## Full-stability, competing-escape, and localized-parent test under one microscopic process

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`
commit `ad15a909`.

Microscopic contract used throughout this report:
the exact leave-one-out reversible finite-N Gibbs heat-bath of report 61.
No empirical-refresh convention is used for the process claims below.

Status:
- pure-k6 scalar branch and its local barrier are exact;
- stationary points / saddle connections reported here are high-precision numerical
  in the full retained X7 space;
- no global capacity theorem is claimed yet.

## 1. Correct full-stability window of the binary k6 prototype

The supercritical parity pitchfork is born at

    g6 = 12/lambda6
       ≈ 5.123427551398616.

The ±k6 daughters are full-X7 minima only until the first k3c transverse
crossing

    g_(3|6)
      ≈ 5.180490619637444.

Therefore the only candidate binary metastability window is

    g6 < g < g_(3|6).

## 2. The k6 minima coexist with a much deeper localized orbit

Full-X7 stationary scans at

    g = 5.13, 5.15, 5.17, 5.179

find two index-zero D12 orbit types:
- the shallow ±k6 orbit;
- the much deeper localized orbit.

Representative potentials:

    g=5.13:
      Phi_k6 ≈ -1.23233e-6
      Phi_loc ≈ -0.795983

    g=5.15:
      Phi_k6 ≈ -2.00497e-5
      Phi_loc ≈ -0.808246

    g=5.17:
      Phi_k6 ≈ -6.13030e-5
      Phi_loc ≈ -0.820521.

Thus local stability of ±k6 does not imply global or even metastable autonomy.

## 3. Competing index-one escape channel

There is an index-one saddle branch connecting the k6 basin to the localized
basin.

Its potential crosses the uniform-saddle level Phi=0 at

    boxed:
    g_iso ≈ 5.150374750802438.

For the k6 minimum let

    B_switch = Phi_uniform - Phi_k6 = -Phi_k6,

and for the observed escape saddle let

    B_escape = Phi_escape_saddle - Phi_k6.

Then:

### for g < g_iso

    Phi_escape_saddle > 0

so

    B_escape > B_switch.

The direct +k6 -> uniform -> -k6 channel is exponentially cheaper than escape
through the observed localized channel.

### for g > g_iso

    Phi_escape_saddle < 0

so

    B_escape < B_switch.

Escape toward the localized phase becomes cheaper than parity switching.

Hence the useful binary window is narrower than the full stability window:

    boxed:
    5.12342755... < g < 5.15037475...

subject to the caveat that a still cheaper undiscovered saddle would narrow it
further.

## 4. Best saddle-level isolation point

Within the two observed competing channels, maximize

    min[
      B_switch,
      B_escape-B_switch
    ].

The optimum occurs at

    boxed:
    g_opt ≈ 5.145228719489142.

There:

    J6 ≈ 0.113032939328268,

    Phi_k6
      ≈ -1.3510972866863e-5,

    Phi_escape_saddle
      ≈ +1.3510972866951e-5,

so

    B_switch
      ≈ 1.3510972866863e-5,

    B_escape
      ≈ 2.7021945733815e-5,

and

    B_escape-B_switch
      ≈ 1.3510972866951e-5.

This is a very small barrier scale.

Ignoring unknown subexponential prefactors, the copy number required for an
exponential margin M is approximately

    N > M / 1.3511e-5.

Thus:

    M=1   -> N ~ 7.4e4
    M=3   -> N ~ 2.22e5
    M=5   -> N ~ 3.70e5
    M=10  -> N ~ 7.40e5.

So a sharply autonomous parity bit appears only at very large N in this
saddle-level estimate.

## 5. Intrawell relaxation

At g_opt the smallest full-X7 Hessian eigenvalue of the k6 minimum gives the
mean-field relaxation rate

    r_slow ≈ 0.00846711747

for unit microscopic refresh clock.

Therefore

    tau_relax ≈ 118.10.

The condition

    tau_relax << tau_switch << tau_escape

requires more than positive barriers:
the switching exponent must also beat log(tau_relax) and the unknown
prefactors.

This makes the large-N requirement stronger, not weaker.

## 6. Direct recursive-child test on the main localized minimum

The main localized minimum has a Z2 reflection stabilizer.

For a representative orientation, the full X7 space splits into

    even dimension 4,
    odd dimension 3.

The localized branch was numerically continued from

    g=3.516

just above the known localization fold, through g=50.

Result:
- the full Hessian remains positive;
- the complete three-dimensional odd block remains positive;
- no reflection-odd soft mode was found.

At g=5.15:

    min odd eigenvalue ≈ 0.1869305891.

At g=3.516:

    min odd eigenvalue ≈ 0.1819046483.

At g=50 the odd eigenvalues approach the positive 1/g scale.

The only soft mode near the birth fold is even, as expected for an ordinary
fold.

## 7. Verdict for the original recursive-pitchfork route

The simple hypothesis

    "each stable localized parent recursively splits into two stable children"

fails for the representative main localized branch over the tested range.

The ±k6 pitchfork remains a genuine endogenous binary-incidence prototype, but:
- it is not local over the main localized phase;
- its autonomy window is narrower than its full stability window;
- useful time-scale separation requires very large N.

Therefore report 131 closes the unrestricted recursive-pitchfork search as a
P0 route.

The next valid route is dynamic coarse-graining of actual metastable basins,
with capacities and memory decay explicitly controlled.
