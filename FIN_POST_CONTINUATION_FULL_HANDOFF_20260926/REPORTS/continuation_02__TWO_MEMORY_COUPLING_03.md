# TWO-MEMORY-COUPLING-03

Status: LIMITED COUNTEREXAMPLE + FIRST-ORDER IDENTIFIABILITY REPAIR.

The product-model separation in TWO-MEMORY-02 is not generically identifiable
under arbitrary weak coupling.

## Counterexample channel

At the uniform visible state let the tree variable y weakly populate the hidden
k=2 simplex direction:

p_i = 1/12 + epsilon y h_{2,psi}(i),
h_{2,psi}(i)=cos(2 theta_i + psi), theta_i=2 pi i/12.

This perturbation leaves all seven retained means equal to zero.  Hence the
heat-bath target remains q=u at that instant, but finite-N jump statistics know
about the full p.

For a rotated k=5 probe

phi_(5,varphi)(i)=sqrt(lambda5/6) cos(5 theta_i+varphi),
F=(sqrt(N) phi^T p)^4,

the exact one-generator excess over the ME7 representative at mu=0 is

(L_full-L_ME)F
 = epsilon y/N * (lambda5^2/3) cos(psi+2 varphi).

For the preregistered pure k=3,4,6 self-quartics the corresponding coefficient
is exactly zero.  Therefore a tree-to-hidden-k2 coupling can imitate the same
"k5 nonzero / k3,k4,k6 zero" selection pattern.  Those controls alone do not
identify the mechanism.

## First-order repair by a preregistered k=5 phase scan

The intrinsic KURTOSIS-THEOREM-11 coefficient is

C_int = lambda5^2(g lambda5-6)/144,

and depends on ||P_H(phi_(5,varphi)^2)||^2, which is independent of varphi.
Thus the intrinsic contribution is constant in probe angle.

The weak external hidden-k2 coupling is a pure second harmonic in varphi.
Pre-register four probe angles

varphi = 0, pi/4, pi/2, 3pi/4.

Their equal-weight average cancels the O(epsilon) coupling exactly, while the
second Fourier harmonic extracts the coupling amplitude and orientation.
Therefore the constant k5 component remains identifiable to first order under
this declared coupling class; residual contamination begins at O(epsilon^2).

If the tree state y is itself resolved from its linear poles/residues, the
cross-spectrum between y and the second-harmonic quartic component supplies an
additional coupling diagnostic.

## Boundary

This is not a generic identifiability theorem for arbitrary nonlinear coupling.
A coupling with an angle-independent O(epsilon) quartic contribution, or an
unobserved mechanism sharing the same invariant component, remains confounded.
