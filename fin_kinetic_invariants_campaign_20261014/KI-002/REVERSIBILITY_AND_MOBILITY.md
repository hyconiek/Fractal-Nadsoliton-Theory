# KI-002: reversible rates, positive mobility, and fluctuation clocks

## Precise finite N generators
For counts n with sum N and any i != j, let m=n-e_i and h_j=(g/N)(A m)_j. The microscopic rate n_i r_ij is defined by

* HB: r_ij=exp(h_j)/sum_k exp(h_k)
* MET: r_ij=(1/12)min(1,exp(h_j-h_i))
* BAR: r_ij=(2/12)[1+exp(-(h_j-h_i))]^-1

These choices are the **K07 g=0 matched distinct-label conventions**. In particular the BAR prefactor 2 is essential. Some earlier source helpers used a different BAR convention: that helper is not the primary input here. All rates at g=0 for any i != j are 1/12, with one nonconstant single-copy Fourier autocorrelation e^-t in the limiting independent-refresh system. Attempts to replace this convention with the N2-sector-normalized clock are unlicensed.

The Gibbs count law is pi_N(n) proportional to [N!/prod n_i!]exp[(g/(2N)) n^T A n]. For n'=n-e_i+e_j, the ratio of count Gibbs probabilities is

pi_N(n')/pi_N(n)=(n_i/(n_j+1)) exp(h_j-h_i)

since A is symmetric with a *constant diagonal*. The common remainder m is the same for the reversed move. HB's forward/reverse conditional-softmax ratio, MET's min/exponential ratio, and BAR's logistic ratio each equal exp(h_j-h_i). Consequently pi_N(n)n_i r_ij(n)=pi_N(n')(n_j+1)r_ji(n') for all three rules. This proves detailed balance of the actual finite-N count chains, not only their continuum approximations.

## Interior force and nonlinear Onsager identity
At N -> infinity the ordered rates are c_ij(p)=p_i rho_ij(p), with h_j=g(Ap)_j. Their two-way ratio obeys

log(c_ij/c_ji)= log(p_i/p_j) + g[(Ap)_j-(Ap)_i] = partial_i F_g - partial_j F_g.

For edge e_ij=e_j-e_i and logarithmic mean Lambda(u,v)=(u-v)/(log u-log v), continuously extended at u=v, set

K(p,g) = sum_{i<j} Lambda(c_ij,c_ji) e_ij e_ij^T.

Then the actual limiting drift is **exactly** b(p,g)=-K(p,g) grad F_g(p) on the tangent. The rate positivity ensures K has kernel span(1) and is strictly positive on the 11-dimensional tangent, for every interior p. Under MET, K is locally Lipschitz, not necessarily C^1, across equal-logit hypersurfaces. The exact identity and local Lipschitz estimates do not require falsely declaring MET C^3.

## Equilibrium linearization and noise
At p_f the constrained free-energy gradient vanishes and c_ij=c_ji>0. Thus

K_f=sum_{i<j} c_ij(p_f)e_ij e_ij^T,
J_f=-K_f H_F on T,
D_f=sum_{i!=j} c_ij(p_f)e_ij e_ij^T=2K_f.

The tangent Hessian H_F has one null and ten positive directions by KI-001. The similarity of -K_f H_F to -K_f^(1/2)H_F K_f^(1/2) on T proves one zero and ten strictly negative rates, with no extra soft mode. The *positive* full tangent K_f is essential to this argument.

Normalize the kinetic left soft vector by

m_K=[r^T K_f,T^(-1)r]^-1,
ell_K^T=m_K r^T K_f,T^(-1),
ell_K^T r=1.

Every inverse here is the true inverse restricted to T (or the corresponding Moore-Penrose inverse in gauge 1^T z=0). Since ell_K^T K_f=m_K r^T on T,

ell_K^T J=0,     ell_K^T D_f ell_K=2m_K.

The leading projected landscape drift is -m_K(a(g-g_f)+(b/2)x_K^2) where x_K=ell_K^T(p-p_f). No independent diffusion coefficient, damping or kinetic timescale has been fitted.

## Heat-bath exact recovery
HB has rho_ij(p_f)=p_f,j, so c_ij=p_f,i p_f,j and K_f=Sigma_f. The inverse-on-T identity gives Sigma_f,T^(-1)r=(Xv modulo a constant), r^T Sigma_f,T^(-1)r=v^T M_f v=1/g_f. Therefore m_HB=g_f, ell_HB=g_f Xv modulo a constant covector, and the admitted SF law is recovered without changing the clock.

## Arithmetic scope
`DERIVED_SOFT_CLOCKS.json` reports the independently cross-checked *nominal* FIN values m_HB≈3.51564471684, m_MET≈0.46077117240, m_BAR≈0.86940912051. No outward interval enclosures for MET/BAR mobility, exact numerical N0, or physical attempt-rate calibration are claimed.
