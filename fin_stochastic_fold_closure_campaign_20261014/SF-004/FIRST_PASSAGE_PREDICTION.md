# SF-004 — Complete canonical first-passage law in a bounded interval

**Status:** analytic prediction conditional on the SF-003 proved *compact-stopped* weak limit, with an independent exact-rational numerical enclosure of one exit-side probability. This does not assume finite-N quantitative TV convergence, global escape or physical clock calibration.

## Frozen observable / no fitted constants

For the FIN fold, the source-oriented coefficients satisfy \(g_f>0\), \(a<0\), \(b>0\). For
\[dZ=-g_f(a\lambda+bZ^2/2)ds+\sqrt{2g_f}\,dB_s\]
put
\[c=(4/b)^{1/3}>0,\quad d=c^2/(2g_f)>0,\quad Y=Z/c,\quad t=s/d,\quad \kappa=-ac/2>0.\]
Brownian scaling gives **exactly**
\[dY_t=(\kappa\lambda-Y_t^2)dt+dW_t.\tag{1}\]
The coefficients 1, -1 and \(\kappa\) follow algebraically; there is no noise amplitude or decay pole chosen to fit a target. This contract freezes \(\lambda=0\), \(Y(0)=0\), canonical interval \((-1,1)\), and canonical horizon \(t=1\). In the original FIN soft clock/coordinate, the same test is \(Z(0)=0\), \(Z\in(-c,c)\), \(s=d\); source fold intervals give roughly \(c\simeq3.227445\), \(d\simeq1.481435\), but exact values are required in the asymptotic theorem.

Let \(\tau=\inf\{t\ge0:|Y_t|\ge1\}\), with \(Y_\tau\in\{-1,+1\}\). The **whole bounded first-passage law** is described by the survival \(u(t,y)=P_y(\tau>t)\) and exit-side finite-time probabilities \(h_\pm(t,y)=P_y(\tau\le t,Y_\tau=\pm1)\); they are nonnegative and sum to 1. In the interior they solve the unique bounded backward parabolic boundary-value problems with generator
\[Gf=\tfrac12f''-y^2f',\]
\[\partial_t u=Gu,\ u(0,y)=1,\ u(t,\pm1)=0,\]
\[\partial_t h_+=Gh_+,\ h_+(0,y)=0,\ h_+(t,-1)=0,\ h_+(t,+1)=1,\]
\[\partial_t h_-=Gh_-,\ h_-(0,y)=0,\ h_-(t,-1)=1,\ h_-(t,+1)=0.\tag{2}\]
The incompatible values at the two corner points \((0,\pm1)\) are interpreted as the usual interior-initial/absorbing-boundary solution; no artificial instant reset or postselection is introduced. Equivalently, the Laplace transforms in time are characterized by **ordinary differential equations**: for each \(q>0\),
\[f_q^\pm(y)=E_y[e^{-q\tau}1_{Y_\tau=\pm1}],\qquad (G-q)f_q^\pm=0\text{ on }(-1,1),\]
with \(f_q^+(-1)=0,f_q^+(+1)=1\) and the reversed conditions for \(f_q^-\); then
\[\int_0^\infty e^{-qt}h_\pm(t,y)dt=\frac{f_q^\pm(y)}{q},\quad \int_0^\infty e^{-qt}u(t,y)dt=\frac{1-f_q^+(y)-f_q^-(y)}q.\tag{3}\]
Equations (2)-(3) uniquely determine **both exit sides and finite-time survival at t=1**; they are exact analytic PDE/ODE representations, not floating numerical predictions. The finite-time probabilities have **not** been approximated here to validated digits; a future table requires a separate deterministic enclosure of the parabolic solver.

### Explicit source-derived exit-side number (at eventual exit, not at t=1)

The eventual upper-side probability \(U(y)=P_y(Y_\tau=+1)\) solves \(GU=0\), \(U(-1)=0,U(+1)=1\). Its derivative satisfies \(U'(y)\propto\exp((2/3)y^3)\). Thus the **exact** answer for the frozen initial condition is
\[U(0)=\frac{\int_{-1}^0 e^{(2/3)y^3}dy}{\int_{-1}^1 e^{(2/3)y^3}dy}.\tag{4}\]
The certificate EXIT_SIDE_CERTIFICATE.json evaluates the integrals by Taylor expansion in \(y^3\), integrating each rational term exactly. With order 25, the absolute remainder for each half-interval is at most \(3(2/3)^{26}/26!\): this follows from \(e^{2/3}<3\) and Taylor's bound. Exact-fraction bounds, followed by outward decimal formatting, give
\[\boxed{0.416855651698382283\ \le U(0)\le\ 0.416855651698382284}.\]
The lower exit probability is its complement, bracketed by \([0.583144348301617716,0.583144348301617717]\). The ratio and all Taylor bounds are free of the FIN fold's decimal midpoint because canonicalization has canceled \(a,b,g_f\) for the specified \(\lambda=0\), canonical boundary/initial conditions.

**Scope warning:** This number refers to the **eventual exit side of a stopped compact diffusion**, not to the finite-horizon probability of exit by canonical t=1, the global physical FIN escape probability, or observed experimental frequency. The finite-horizon law is specified by the exact (2)-(3) representation; values would need a separately certified BVP solve to print as digits. Weak convergence transfers finite-horizon exit event probabilities (and their side) via the boundary continuity proved in SF-003; it does not supply a finite-N speed or convergence of unbounded physical residence-time moments.
