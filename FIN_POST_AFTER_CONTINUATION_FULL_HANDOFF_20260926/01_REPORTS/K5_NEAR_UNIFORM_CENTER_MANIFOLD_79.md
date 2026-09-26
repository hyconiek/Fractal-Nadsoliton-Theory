# K5-NEAR-UNIFORM-CENTER-MANIFOLD-79
## Degree-12 angular anisotropy explains the numerically delicate k=5 branch cluster

Date: 2026-09-26

Status:
- exact symmetry order and exact positive quartic coefficient;
- high-precision numerical center-manifold estimate of the first angular coefficient;
- no interval certificate for the degree-12 coefficient yet.

## 1. Critical sector

At the uniform state the k=5 pair becomes soft at

    g5 = 12/lambda5
       ≈ 5.220554796949195.

The critical coordinate is a faithful real two-dimensional D12 representation,

    z = r exp(i phi).

Because 5 is coprime to 12, this representation sees the full twelvefold
rotation structure.

## 2. Allowed invariants

D12 invariance permits the radial invariant

    r^2

and radial powers thereof.

The first nontrivial angular invariant is

    Re(z^12)
      = r^12 cos(12 phi).

Therefore every angular derivative of the reduced potential through degree 11
vanishes identically.

This is why the k=5 sector is numerically almost O(2)-symmetric close to the
uniform threshold.

## 3. Quartic radial coefficient

At uniformity, for a unit k=5 Fourier direction f,

    D^4 Phi[f,f,f,f]
      = -kappa4(f)
      ≈ 0.0550374041046
      >0.

Thus the coefficient of r^4 in the reduced potential is

    a4 = D^4 Phi / 24
       ≈ 0.00229322517102
       >0.

So the primary radial instability is locally stabilizing/supercritical.

The finite-amplitude branches already present at g=5 are therefore not caused
by a negative local quartic at uniformity.

## 4. Stable-mode slaving

At g=g5, fix the k=5 center coordinate `(r,phi)` and solve the five complement
stationarity equations for the retained modes

    k=3, k=4, k=6.

The first forced complement orders are consistent with Fourier arithmetic:

    k=3 appears at O(r^3),
    k=4 appears at O(r^4),
    k=6 appears at O(r^6).

After this slaving, compare the two inequivalent reflection axes

    phi=0
and
    phi=pi/12.

The energy difference has the form

    Phi_eff(r,0)-Phi_eff(r,pi/12)
      =
      2 C12 r^12
      +O(r^14).

High-precision solves for r down to 0.002 and extrapolation in r^2 give

    boxed:
    C12 ≈ 0.00252080442115
    >0.

## 5. Which reflection axis is lower?

Because C12>0,

    cos(12phi)=-1

is the lower-energy angular family near uniformity.

Thus the preferred reflection axes are

    phi = (2n+1)pi/12,

while

    phi = n pi/6

are the higher angular stationary axes.

The two sets are inequivalent D12 orbits, each with reflection stabilizer Z2.

## 6. Angular stiffness scale

For

    Phi_aniso=C12 r^12 cos(12phi),

the tangential Hessian eigenvalue scales as

    lambda_ang
      ~ 144 C12 r^10

on a stable angular axis.

Thus

    boxed:
    lambda_ang = O(r^10).

This extremely high power explains why:
- continuation near g5 is ill-conditioned;
- apparent angular zero modes persist over a visible parameter interval;
- tiny numerical perturbations can switch between nearby symmetry-related
  branches.

The small eigenvalues are structural, not evidence of a new continuous
symmetry.

## 7. Identification of the known reflection branch

The reflection-symmetric component containing atlas orbits i=4 and i=13 has
k=5 phase satisfying

    cos(12phi)=-1.

Its small-amplitude sheet therefore lies on the lower-energy primary k=5 axis
predicted by the center-manifold coefficient.

Numerical continuation drives that sheet toward the uniform threshold
`g5=12/lambda5`.

## 8. Secondary structure

Farther from uniformity, higher-order radial/angular terms can reverse angular
stability.

Indeed the same reflection branch develops additional odd-sector zero
eigenvalues at finite amplitude. The first such event relevant to the g=5
atlas is analyzed in report 80.

## 9. Boundary

C12 is currently a high-precision numerical center-manifold coefficient, not
an interval-certified constant.

The exact conclusions are:
- angular anisotropy cannot occur before degree 12;
- the quartic radial coefficient is positive.

The numerical conclusion is the sign/magnitude of the fully slaved degree-12
coefficient.
