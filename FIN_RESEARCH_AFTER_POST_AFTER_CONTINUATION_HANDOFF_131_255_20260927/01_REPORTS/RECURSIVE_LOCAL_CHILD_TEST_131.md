# RECURSIVE-LOCAL-CHILD-TEST-131
## A representative localized FIN minimum does not reproduce the k6 binary pitchfork mechanism

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`
commit `ad15a9098ecc5e1282f964ea8b159a8ec608d7c5`.

Status:
- microscopic convention frozen to the exact leave-one-out Gibbs heat bath;
- exact symmetry criterion for the candidate Z2 child mode;
- dense numerical continuation of the localized stable branch;
- high-g sign supported by the one-label concentration asymptotics;
- NOT yet an interval cover of the whole finite-g branch.

## 1. Microscopic contract

All metastability and process-level statements in this campaign use the exact
leave-one-out Gibbs heat bath of report 61:

    q_j^(-i)(n)
      =
      softmax_j[
        (g/N) A7 (n-e_i)
      ],

with count jump

    n -> n-e_i+e_j

at rate

    n_i q_j^(-i)(n).

Its stationary count law is the multinomial degeneracy times the pure
finite-copy Gibbs weight.

The older empirical-refresh convention is retained only as a comparison model.

## 2. Convention transfer rule

The exact leave-one-out drift has expansion

    b_N(p)
      =
      q(p)-p
      -(g/N) S_q A7 p
      +O(N^-2),

where

    q=softmax(g A7 p),
    S_q=diag(q)-q q^T.

The leading jump covariance is

    D_N(p)
      =
      (1/N)
      [
        diag(p)+diag(q)-p q^T-q p^T
      ]
      +O(N^-2).

At equilibrium p=q=p*,

    D_N(p*)
      =
      (2/N) S(p*)
      +O(N^-2).

Therefore:
- the deterministic N->infinity heat-bath lane is common;
- the leading Gaussian/CLT mobility-noise structure is common;
- O(N^-1/2) hidden-preparation effects built only from that leading Gaussian
  lane survive the convention change at leading order;
- O(1/N) Edgeworth, stationary-bias, and finite-N invariant-measure claims
  must be recomputed for leave-one-out and must not be copied from empirical
  refresh.

An independent numerical replay gives N^2 times the first-order drift remainder
approximately constant:
    N=20   -> 0.04336
    N=40   -> 0.04319
    N=80   -> 0.04309
    N=160  -> 0.04303.

## 3. Localized parent and its only canonical binary symmetry channel

Take the main localized index-zero branch born at the certified simple fold

    g_f≈3.51564471684.

A representative can be chosen reflection-even with only real
k3,k4,k5,k6 coordinates.

Its residual generic stabilizer is Z2 reflection.

A local canonical binary child pair obtained by breaking that residual
reflection would require a zero eigenvalue in the reflection-odd block

    H_odd = H|_(k3s,k4s,k5s).

So a necessary local pitchfork condition is

    lambda_min(H_odd)=0.

## 4. Continuation result

The stable localized branch was continued from just above the lower fold
through g=1000.

Representative values:

    near fold, g=3.51565:
      lambda_min(H_odd)
        ≈ 0.1798878404;

      g lambda_min(H_odd)
        ≈ 0.6324226861.

    g=5:
      lambda_min(H_odd)
        ≈ 0.1913061655.

    g=1000:
      lambda_min(H_odd)
        ≈ 0.001,
      so
      g lambda_min(H_odd)
        ≈ 1.

Across the numerical continuation no odd eigenvalue approaches zero.

The only zero encountered at the lower endpoint is the known EVEN fold mode.

At large g the branch converges to a one-label support.  Its categorical
covariance vanishes exponentially while the dual Hessian tends to

    H ~ I/g,

so the odd block is asymptotically positive.

## 5. Verdict for the simplest recursive-pitchfork hypothesis

There is no evidence that the representative localized parent reproduces the
uniform -> ±k6 mechanism by locally breaking its residual Z2 symmetry.

The bounded test therefore returns

    boxed:
    FAIL_NO_LOCAL_Z2_CHILD_ON_LOCALIZED_BRANCH

for the simplest recursive binary-refinement route.

This does NOT prove that no other parameter, nonlocal bifurcation, or
non-binary metastable decomposition exists.

It does mean that the next campaign should not continue by blindly searching
for another copy of the k6 pitchfork on this branch.

## 6. Proof obligation if a theorem is desired

To upgrade the negative result to a full theorem:
- interval-cover the reflection-odd block on a finite branch segment
  [g_f+epsilon,G];
- combine it with an analytic high-g covariance bound for g>=G.

That is much smaller than another global stationary atlas.
