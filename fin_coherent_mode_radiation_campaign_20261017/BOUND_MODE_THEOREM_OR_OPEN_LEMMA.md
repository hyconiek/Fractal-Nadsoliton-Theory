# RP-002 — isolated transverse mode in one frozen source-compatible coefficient family

**Disposition: PROVED_CONDITIONAL_BOUND_MODE_ON_DECLARED_FAMILY**, not a proof about naturally selected FIN parameters. This file is a mathematical argument, not independent peer review or formal proof verification.

## Frozen family (before spectral and on-shell evaluation)

Use exactly the existing *form* of the 1+1 canonical radial action: one real phi, 12 complex Psi components, real symmetric 12x12 K_s with zero diagonal, diagonal M2 and Y, and arbitrary finite **unchanged signed** Psi quartic and sextic self coefficients. No phi^6 term. In dimensionless model coordinates choose the phi4 member v=k=1, U(phi)=(phi^2-1)^2/2, K(x)=tanh x, L_phi=-d_x^2+4-6 sech^2 x.

The parameter family F_RP is specified **once** as follows (mathematical witness, not a FIN fit):

* a in [-3/4,-1/4]; M2=diag(a,a+Delta_2,...,a+Delta_12), each Delta_j in [0.10,0.20]; all diagonal masses are strictly negative.
* Y=diag(3,3+d_2,...,3+d_12), |d_j|<=0.01. Require d_2 != d_3.
* K_s has only K_23=K_32=h nonzero, 0<|h|<=0.02. Hence [M2+K_s,Y] !=0 when d_2 != d_3. First complex channel e1 is uncoupled **at the linear Hessian only**; the other 11 complex channels and all velocities remain in the PDE. K_s has the source-bound zero diagonal.
* All remaining already sourced radial Psi quartic/sextic coefficients are carried symbolically (finite, fixed, no new sign selections). Their Hessian at Psi=0 vanishes, and for sufficiently small full data NS's local polynomial energy argument applies. One is **not** asserting the actual source chooses F_RP.

Put A=M2+K_s and B=2Y. On each 12-real copy L_perp=-d_x^2 I + A+B tanh^2 x; the second real copy is identical. The first component decouples from the remaining real 11 component matrix but there is **no physically licensed real-only truncation**; f and i f form the explicit complex doublet.

## Exact eigenfunction and gap

Let e1 be the first coordinate vector. Define f(x)=sqrt(3/4) sech^2 x e1. The identity int_R sech^4 x dx=4/3 proves ||f||_2=1. The Poeschel–Teller factorization with s=2 gives

    L_perp f=(a+2) f = omega^2 f,
    omega=sqrt(a+2), omega^2 in [5/4,7/4].

f is even, real, positive in e1, and |f(x)| <= 4 sqrt(3/4) exp(-2|x|). There is an exactly two-dimensional real phase eigenspace {f, i f}; the real 12-block's ground mode is simple. The e1 excited PT eigenvalue is a+5, and its threshold is a+6.

For the complement use Y_comp >= 2.99 I, so B_comp>=5.98 I and M2_comp+K_s >= (a+0.10-0.02) I=(a+0.08)I, regardless of commutation. Thus

    L_comp >= [-d_x^2+a+0.08+5.98 tanh^2 x] I.

The scalar lower comparison has lowest eigenvalue a+0.08+s_-, s_-=(sqrt(1+4*5.98)-1)/2 >1.995. Hence all complement spectral values are >a+2.075; in particular none coincides with a+2. This pays real-block simplicity without omitting the complement. The asymptotic complement vacuum matrix is >=(a+0.10+5.98-0.02)I=a+6.06, while e1 has vacuum threshold a+6. Standard compactness of the decaying matrix perturbation implies sigma_ess(L_perp)=[a+6,infinity). Thus a+2 is isolated below the full essential edge by 4.

NS positivity for all real and imaginary blocks: V_inf>=a+6 I and b=max(0,lambda_max B)<=6.02. Therefore the imported NS quantitative bound is delta_Psi>=a+6-s_+^2, s_+=(sqrt(1+4*6.02)-1)/2<2.005, so delta_Psi>1.22; scalar gap modulo K' is exactly delta_phi=3. Consequently NS's full field small-data global orbital result applies conditionally under its Cauchy/continuation assumptions, with signed higher Psi interactions unchanged. This gives persistence, *not* an outgoing flux.

## Scope

The family permits a source-type nonzero, noncommuting K_s in channels 2 and 3, but uses a decoupled first **linear** channel for an exact solvable witness. It is an open mathematical family of freely supplied coefficients, not a provenance certificate for actual FIN coefficients. Other parameter tuples may lack a localized mode; a positive Hessian gap alone never establishes a mode.
