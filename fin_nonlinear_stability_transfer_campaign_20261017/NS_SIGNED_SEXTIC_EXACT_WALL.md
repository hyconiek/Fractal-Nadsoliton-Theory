# NS-001/002 exact signed-sextic wall and a fully algebraic subregion

This is a NEW conditional analytical calculation in the explicitly declared standard even sextic U; it is *not* a claim that this scalar formula exhausts the inaccessible FIN action. It neither assigns nor fits the existing coefficients. If the supplied legacy radial U has a different normalization or extra mixed powers, substitute its actual formula before applying the certificate.

## Keep the signed mass convention (no benign replacement)

Take U(phi)=-m^2 phi^2/2 + lambda phi^4/4 + g phi^6/6 + const with lambda>0, v>0, and g *signed, possibly nonzero*. The vacuum stationarity constraint is m^2=lambda*v^2+g*v^4, with the constant chosen to make U(+-v)=0. This gives the **exact factorization**

U(phi)=(phi^2-v^2)^2 [lambda/4 + (g/6)(phi^2+2v^2)].

Write C=lambda/2+2g v^2/3, D=g v^2/3, k=v sqrt(C+D)=v sqrt(lambda/2+g v^2), theta=D/C. If C>0, C+D>0 and |theta|<1, the explicitly stationary monotone kink is

K_g(x)=v sqrt(C) tanh(kx) / sqrt(C+D-D tanh^2(kx))
      =v tanh(kx) /sqrt(1+theta sech^2(kx)).

Indeed K_g(+-infty)=+-v, K_g'>0 and the exact first integral is K_g'^2 =2 U(K_g). This is the full sextic kink with its *nonzero g*; a bare tanh is not exact unless g=0. No attempt is made to identify its width k with SI length.

## Fully checkable error bounds (no new potential)

Put t=tanh(kx), q=1-t^2. Then

K_g^2/v^2=t^2/(1+theta q),
|K_g^2/v^2-t^2|=|theta| t^2 q/|1+theta q| <= |theta|/[4(1-|theta|)]=r_theta.

Thus the true full transverse operator differs pointwise from the tanh comparison by at most ||B||_op r_theta, B=2v^2Y. No commutativity is required.

The scalar Hessian deviation from H0=-d^2+4k^2-6k^2 sech^2(kx) is exactly

U''(K_g)-(4k^2-6k^2 sech^2(kx))
  =g v^4 [1-6t^2+5 (K_g/v)^4]
   +3 lambda v^2 [(K_g/v)^2-t^2].

Since |K_g/v|<=1, a completely explicit (conservative) sufficient error bound is

e_6 :=12 |g| v^4 +3 lambda v^2 r_theta.

Define mu=lambda_min(A+B), b=max(0,lambda_max(B)), s=(sqrt(1+4b/k^2)-1)/2. Then **for every tuple** in the region

P_6 = {v,lambda>0, g real, m^2=lambda v^2+g v^4,
       C>0,C+D>0, |theta|<1,
       delta_phi :=3k^2-e_6>0,
       delta_Psi :=mu-k^2 s^2-||B||_op*r_theta>0,
       all coefficients/signatures exactly as frozen in NS_CONTRACT},

all scalar and 24 real transverse coercivity and the nonlinear orbital theorem of NS_ORBITAL_THEOREM hold (subject to the exact unchanged R4 polynomial and Cauchy hypotheses). The inequalities are checkable using only the full frozen coefficients, not hand-picked target-particle masses.

## Nonemptiness with g != 0 and a signed negative core mass

Choose any mathematical control with lambda,v>0 and an indefinite/non-positive A, noncommuting B, for which the g=0 PT comparison has a strict positive margin. This is possible for the *inequality class*: with fixed nonconstant positive D and fixed symmetric E satisfying [D,E]!=0, A=-D, B=t I+E has mu-k^2s^2>0 for some sufficiently large but **finite** t. Once such a finite tuple is frozen, every inequality defining P_6 depends continuously on g near zero. Therefore there is a nonempty two-sided interval of **nonzero signed g** sufficiently close to zero in which P_6 remains satisfied; it is an open region with noncommuting matrices and negative A at the wall centre. This existential control neither changes a FIN coefficient after observing data nor proves that its fixed K_s/Y source allows this region.

The actual FIN-specific result remains OPEN_PARAMETER_MEMBERSHIP until its *actual* signed scalar and 12x12 arrays are evaluated. A failed sufficient region is not a discovered instability. The existence theorem is for a specified predeclared mathematical family, not the observed universe.
