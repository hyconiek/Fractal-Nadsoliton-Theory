# SRN-003 — Analysis-only mollifier and discrete OU corrector

In observed coordinates w=(w1,w3) define H(w)=1{|w3|>|w1|}; let rho in C_c^infty(R²), rho>=0, integral rho=1, support inside unit ball; rho_eta(w)=eta^-2 rho(w/eta), and H_eta=H*rho_eta. H is constant on each connected component off the two diagonal boundary lines; hence |H_eta-H|<=1_{dist(w,L_+ union L_-)<=eta} up to a fixed universal dilation constant. For multiindices beta with |beta|=j<=3,
 ||D^beta H_eta||_infty<=C_j eta^-j.
No random noise is applied to raw FIN; H_eta is never seen by the classifier.

For r in {1,3}, (Fr,Rr) follows the **analysis** OU generator with gamma1=1, gamma3 in the two fixed positive values and private diffusion coefficient (1-alpha), 0<alpha<=1/2. The observed linear map is w_r=sqrt(alpha)Fr+Rr, whose norm remains bounded uniformly in alpha. Let P_t be this 4D OU semigroup, and mu_eta=E_G H_eta(W1,W3). The corrector
 psi_eta(x)=sum_{l>=0} (P_l(H_eta o L_alpha)(x)-mu_eta)
converges locally uniformly; (I-P_1)psi_eta=H_eta o L_alpha-mu_eta. Coupling an OU started at x to an independent stationary OU gives
 |P_t H_eta(x)-mu_eta|<=C eta^-1 e^-lambda t (1+||x||), lambda=min(1,gamma3)>0.
Hence |psi_eta(x)|<=C eta^-1(1+||x||). By differentiating the explicit linear OU transition,
 ||D^j P_t H_eta||_infty<=C eta^-j e^-j lambda t, j=1,2,3.
Consequently for j=1,2,3, sup_{0<=s<=1}||D^j P_s psi_eta||_infty<=C eta^-j. All constants uniform in alpha<=1/2 and in the two fixed gammas. No global zero-order sup bound is asserted or needed.

Reflected stationary hypergeometric/private and MT fast moments give E_piD ||phi||²<=C, so E_piD psi_eta(phi)²<=C eta^-2. Duhamel applied to the full finite count chain (not an assumed Markov quotient phi) yields
 E_piD |P^D_1 psi_eta(phi)-(P_1 psi_eta)(phi)|<=C eta^-3 epsilon_N,k,
 epsilon_N,k=N^-1/4+k^-1/2+(N/sqrt(k))pi_D(B_D),
for eta<=1. The order eta^-3 pays *all* Taylor derivatives rather than hiding smoothing cost.
