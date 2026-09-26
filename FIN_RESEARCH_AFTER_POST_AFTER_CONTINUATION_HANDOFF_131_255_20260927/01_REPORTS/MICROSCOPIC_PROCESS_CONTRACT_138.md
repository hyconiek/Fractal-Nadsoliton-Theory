# MICROSCOPIC-PROCESS-CONTRACT-138
## One declared microscopic process and a transfer dictionary for older FIN results

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, commit
`ad15a9098ecc5e1282f964ea8b159a8ec608d7c5`.

Primary microscopic process for the current campaign:

    exact leave-one-out Gibbs heat-bath.

The older empirical-refresh process is retained only as an explicitly separate
comparison model.

## 1. Objects shared at leading mean-field level

Both conventions have the same declared large-N deterministic lane:

    V_g(p)
      = D(p||u0)
        -(g/2) p^T A7 p,

    q(p)
      = softmax(g A7 p),

    p_dot
      = q(p)-p

up to the chosen global refresh clock.

Therefore the following are properties of the common declared mean-field model,
not distinguishing finite-N conventions:
- stationary equation p=q(p);
- stationary branch landscape;
- Hessian/soft-mode calculations;
- the exact mean-field Lyapunov identity;
- the associated Onsager representation of that deterministic flow.

## 2. Exact finite-N equilibrium is convention-specific

### leave-one-out

The exact count process has the pure Gibbs occupation weight

    pi_LOO(n)
      proportional to
      (N!/prod_i n_i!)
      exp[(Ng/2) p^T A7 p].

This process is exactly reversible and is the process used in reports 131+.

### empirical refresh

The predecessor empirical-refresh count chain has an additional finite-N factor

    Z_g(p).

Therefore:

    pi_empirical != pi_LOO

at finite N even though their mean-field drift agrees.

Any finite-N statement sensitive to the stationary measure, capacities,
committors, prefactors or O(1/N) coefficients must name the convention.

## 3. What transfers safely

### Safe as common mean-field mathematics
- stationary branches and local bifurcations of V_g;
- large-N rate-function barrier HEIGHTS at O(N), when the extra stationary
  factor is only subexponential in N;
- deterministic Onsager/Lyapunov structure.

### Safe only after a leading-generator check
- Gaussian fluctuation/FDT coefficients;
- leading hidden OU decomposition;
- O(N^-1/2) nonstationary preparation effects.

The two finite-N chains have the same leading mean-field law, but the exact
jump-rate expansion should be checked before treating every fluctuation
coefficient as identical.

## 4. What does NOT transfer automatically

The following must be rederived for leave-one-out before being used in the new
multiscale campaign:
- stationary O(1/N) Edgeworth coefficients previously derived for the
  empirical-refresh invariant;
- finite-N projector corrections;
- non-Gaussian finite-N hidden cumulants;
- exact finite-N memory kernels;
- capacities and Eyring-Kramers/potential-theory prefactors;
- multi-time coarse-grained process laws.

Reports 131-143 therefore compute their finite-N process quantities directly
from the leave-one-out generator.

## 5. Campaign rule

Every new dynamical result must carry one of the tags:

    MEAN-FIELD-COMMON

or

    LEAVE-ONE-OUT-EXACT

or

    EMPIRICAL-REFRESH-ONLY.

No formula is transferred between the last two categories without an explicit
generator-expansion comparison.

This prevents convention drift from masquerading as FIN universality.
