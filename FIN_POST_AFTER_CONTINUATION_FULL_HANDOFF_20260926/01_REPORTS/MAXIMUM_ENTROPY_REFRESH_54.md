# MAXIMUM-ENTROPY-REFRESH-54
## Complete refresh is the unique maximum-conditional-entropy update at fixed target law

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact finite-state information-theoretic theorem;
- supplies a naturality principle for the *shape* of the one-event heat-bath
  update;
- does not determine the target q or the event clock rate.

## 1. Setup

Let q be an interior probability distribution on Q finite labels.

Let K be a row-stochastic one-event transition kernel and suppose q is
stationary:

    q^T K = q^T.

Take X distributed as q and Y generated from K:

    P(X=i,Y=j)=q_i K_ij.

Because q is stationary,

    Y ~ q.

## 2. Conditional entropy bound

The mutual information identity gives

    I(X;Y)=H(Y)-H(Y|X) >= 0.

Since Y~q,

    H(Y)=H(q).

Therefore

    H(Y|X) <= H(q).

Equality holds iff

    I(X;Y)=0,

that is, iff X and Y are independent.

Because q_i>0 for every i, independence implies

    q_i K_ij = q_i q_j

for every i,j, hence uniquely

    boxed:
    K_ij=q_j.

So the complete-refresh kernel

    K = 1 q^T

is the unique maximizer of conditional entropy among q-stationary one-event
kernels.

## 3. Interpretation

The maximizing update keeps no label information from the previous state beyond
the externally supplied target q.

Equivalently:

    old label -> discard -> draw new label independently from q.

Thus complete refresh is not merely a convenient reversible kernel. It is the
unique maximum-forgetting event at fixed target law.

## 4. Continuous-time embedding

Let refresh events occur with Poisson rate rho.

Then

    Q = rho(K-I)
      = rho(1 q^T-I).

All q-centered observables decay with the same rate rho.

The entropy principle selects the generator shape but not rho.

Multiplying rho by a positive constant only changes the clock scale.

## 5. Why the event formulation matters

If one tries to maximize a continuous-time entropy *rate* without fixing an
event-rate budget, the objective is scale-unbounded because all transition
rates can be multiplied by a constant.

The theorem therefore concerns a normalized one-event update. A separate
Poisson/event-rate parameter supplies the time scale.

## 6. Boundary cases

If some q_i=0, rows associated with zero-probability starting states are not
constrained by the q-weighted conditional entropy. Uniqueness is therefore
claimed for interior q, which is exactly the regime used by the accepted
interior heat-bath expansion.

## 7. FIN consequence

Once a FIN rule supplies the target q(mu), maximum conditional entropy selects
the microscopic complete-refresh event uniquely.

Still not supplied:
- why the target is q(mu);
- the coupling g entering that target;
- the physical rate rho;
- physical seconds.

The next task is to test whether q(mu) itself follows from the existing FIN
entropy/variational structure.
