# NS-002 — noncommuting quantitative transverse and scalar spectral gaps

## Exact frozen, open coefficient region P_NS

Given the single action in NS_CONTRACT.md and its **actual** heteroclinic K, choose the comparison inverse length k>0 once, before evaluating either test. Write:

- A=M_2+K_s (real symmetric 12x12), B=2v^2 Y, V_inf=A+B;
- mu=lambda_min(V_inf), b=max(0,lambda_max(B)), s=(sqrt(1+4b/k^2)-1)/2;
- rho(x)=K(x)^2/v^2, e_wall=sup_x |K(x)^2-v^2 tanh^2(kx)|, and e_scalar=sup_x |U''(K(x))-[4k^2-6k^2 sech^2(kx)]|;
- e_trans=(||B||_op/v^2) e_wall, delta_Psi=mu-k^2*s^2-e_trans;
- delta_phi=3k^2-e_scalar.

The *predeclared* region P_NS is all tuples with v,k>0, the exact stationary monotone wall, finite e_wall and e_scalar, delta_Psi>0 and delta_phi>0, and locally C^3 radial V with a well-posed canonical wave equation. This is an explicit **sufficient** region in norms/eigenvalues of frozen data, not an inferred natural parameter value; it is not claimed to contain actual FIN parameters. All 12 complex transverse modes (24 real) are included. Crucially, for non-tanh sextic walls the e_trans term cannot be dropped.

## Proof of noncommuting transverse gap

A+2K^2Y = V_inf-B sech^2(kx) + B*(K^2/v^2-tanh^2(kx)).
The final multiplication operator has operator norm <=e_trans. Since B <= b I and V_inf >= mu I independently of commutation, the associated quadratic form satisfies

Q_Psi[z]>= <z,[-d_x^2+mu-b sech^2(kx)] z> - e_trans ||z||_2^2.

The scalar Poeschel--Teller operator -d_x^2-b sech^2(kx) has ground eigenvalue -k^2 s^2 (s(s+1)=b/k^2). Tensor this scalar form with R^24. Therefore Q_Psi>=delta_Psi ||z||_2^2. This is a uniform *strictly positive* gap for EVERY coefficient tuple in P_NS without simultaneous diagonalization. No finite spectral digit sweep is used.

Alternative geometry-only sufficient test: if |K|<=v and A>=alpha I, V_inf>=beta I with alpha,beta>0, then A+2K^2Y=(1-rho)A+rho V_inf>=min(alpha,beta) I, so Q_Psi>=min(alpha,beta)||z||_2^2. This also needs no commutativity, but may fail for signed tachyonic A; it is not substituted for the main delta_Psi test.

If delta_Psi<=0 or the exact e_trans is uncomputed, this sufficient test is *inconclusive*; it neither proves negative spectrum nor instability. The independently supplied UE Rayleigh negative-direction criterion applies only to an exact tanh wall (or with a controlled wall-error correction).

## Scalar gap and translation direction

Since K' >0 and L_phi K'=0, Sturm oscillation/factorization implies zero is the simple ground eigenvalue of L_phi and there is no negative scalar direction. H0=-d_x^2+4k^2-6k^2 sech^2(kx) has eigenvalues 0 and 3k^2 below continuum 4k^2. The bounded multiplication difference between L_phi and H0 is <=e_scalar. The min-max principle implies its second spectral level >=3k^2-e_scalar. Thus

Q_phi[h]>=delta_phi ||h||_2^2 for h perpendicular to K'.

This step compares to the correct scalar phi4 **control**, without identifying the actual quintic-force sextic wall with tanh. The condition e_scalar<3k^2 is a stringent sufficient test, not a necessary condition.

## H1 coercivity with all modes

Let delta=min(delta_phi,delta_Psi)>0 and C_V=max(0, 2k^2+e_scalar, b-mu+e_trans). Pointwise Hessian lower bounds give Q[h,z] >= ||(h',z')||_2^2-C_V||(h,z)||_2^2. Combined with the spectral L2 gap on the codimension-one subspace <h,K'>=0 yields

Q[h,z] >= c_H (||h||_{H1}^2+||z||_{H1}^2), c_H=delta/(delta+1+C_V)>0.

This is coercivity of the **full 25-real-position** Hessian, with exactly one translation zero mode, and of all 25 associated velocity components through kinetic energy. No imaginary Psi block is silently removed.

## Region nonemptiness versus source

At the phi4 control point e_wall=e_scalar=0. The transverse gap region is mathematically nonempty even when the signed A matrix has negative eigenvalues: for a fixed real diagonal nonconstant positive D, a fixed real symmetric E with [D,E]!=0 and sufficiently small nonzero ||E||, take A=-D, B=tI+E (t>0 treated as a free *mathematical* family parameter). Then mu>=t-||D||-||E||, b<=t+||E|| and delta_Psi>=k sqrt(t)+O(1)>0 as t->infty. It is open under sufficiently small changes in its coefficients. This is a proof of feasibility of the **inequality region only**, not a choice of FIN couplings or particles, and it does not show that actual kernel restrictions permit this family. If U has a nonzero sextic coefficient, nonemptiness of its scalar close-to-phi4 subregion requires showing the exact quadrature wall satisfies the two wall-error bounds. That input is not independently paid for the actual coefficient tuple.

## Explicit nonzero signed phi6 refinement (paid)

For the exact signed sextic polynomial U=-m²phi²/2+lambda phi⁴/4+g phi⁶/6+const with m²=lambda v²+g v⁴, the full kink K_g and *explicit* upper bounds e_wall/v² <= r_theta and e_scalar <= e_6 are derived in NS_SIGNED_SEXTIC_EXACT_WALL.md. Consequently the fully algebraic region P_6 (including small signed NONZERO g) is a genuine nonempty sufficient class of conditional nonlinear stability. This refinement supersedes any use of an assumed tanh kink at nonzero g and does not assign it to unverified FIN coefficient bytes.
