# STRICT-RATIO-SOURCE-52
## Conditional Fisher whitening does not fix the remaining strict kinetic ratio

Date: 2026-09-26

Status:
- exact covariance/OU no-go;
- exact direct-strict-family no-go;
- exact conditional statement showing hidden rate isotropy would force rho=0,
  but that isotropy is an additional kinetic premise.

## 1. Remaining ratio

After the optional acceptance cocycle of report 51, a direct strict-side
generator can be written, up to an overall clock rate, with

    rho=eta/gamma.

The question is whether Fisher/Schur whitening or existing static data fixes
rho.

## 2. Whitening is not kinetic isotropy

Let the conditional hidden covariance be Sigma>0.

Whiten

    q=Sigma^{-1/2}z.

Stationary covariance is now I.

But for ANY symmetric positive-definite matrix C,

    dq = -C q dt + sqrt(2C) dW

has stationary covariance I because the Lyapunov equation is

    C I + I C = 2C.

Therefore whitening fixes the covariance metric but leaves the full SPD
relaxation matrix C free.

So:

    Fisher whitening != equal relaxation rates.

## 3. Unwhitened form

Equivalently, for any SPD mobility M, the Gaussian law N(0,Sigma) is reversible
for

    dz = -M Sigma^{-1} z dt + sqrt(2M)dW.

The stationary covariance is Sigma for all M.

Thus the statistical geometry alone cannot source the kinetic mobility.

## 4. Direct strict family

At uniformity,

    Q=-gamma P_C-eta A

has hidden Fourier rates

    r_1=gamma+eta lambda_1,
    r_2=gamma+eta lambda_2.

Since strict FIN has

    lambda_1 != lambda_2,

requiring equal k=1 and k=2 hidden rates forces

    eta=0.

Hence

    rho=0

would follow IF one separately imposes isotropic hidden decay.

But that requirement is exactly a kinetic condition; it does not follow from
the hidden covariance being whitenable.

## 5. Localized equilibrium also leaves rho free

Using the half-density detailed-balance deformation,

    q_ij =
      (gamma/12+eta W_ij)
      sqrt(p_j/p_i),

every eta>=0 gives:
- the same prescribed localized equilibrium p;
- positivity;
- reversibility;
- D12 covariance;
- the same half-density acceptance rule.

Changing eta changes the relaxation spectrum while leaving the equilibrium
state unchanged.

Thus equilibrium + detailed balance + the optional cocycle do not fix rho.

## 6. Relation to the declared heat-bath model

The admitted finite-N heat-bath refresh chooses no direct `eta A` mobility
term: in this parameterization its microscopic direct ratio is rho=0.

But that is inherited from the supplied refresh law.

It is not derived by Fisher geometry.

## 7. Disposition

    FISHER/SCHUR_GEOMETRY_DOES_NOT_FIX_RHO.

The heat-bath value rho=0 is recovered only after importing its isotropic
refresh dynamics.

The next task is to understand how strict spectral information nevertheless
enters the effective heat-bath rates.
