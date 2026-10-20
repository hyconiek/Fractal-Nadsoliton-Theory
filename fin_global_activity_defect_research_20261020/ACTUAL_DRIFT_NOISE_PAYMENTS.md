# AD-003: the two original-source critical payments are now proved

Date: 2026-10-20. Conditional on the actual inherited static/local inputs.
All expectations are under exactly pi_D,N for the reflected comparison.

Put epsilon=1/N, L_N=1+log(N+1), eta_N=pi_D(B_D),
f(z)=q(chi lambda z-z^3/3), and use the actual beta_N,G_N from the IS
discrete interface. Keep the frozen d_N=1+N(q-1) normalization.

## 1. Actual source versus derived auxiliary off-diagonal flow

Let beta_app,G_app use off-diagonal P_app(x,y)=u_xy/pi_x, restricted
to internal D edges exactly as for the true source. At critical scaling
h_N=1/(d_N sqrt(N)), each Z jump is at most 2N^-3/4. Therefore (4) of
AD-002 with w=1+|Z| gives

    E_piD[(1+|Z|)|beta_N-beta_app|]
       <=2d_N N^-1/4 sum_inside (1+|Z_x|)|c_xy-u_xy|/pi(D)
       <=C N^3/4 sqrt(L_N[N^-2+eta_N]).              (1)

With weight one and squared jumps,

    E_piD|G_N-G_app|
       <=(4d_N/N)sum_inside |c_xy-u_xy|/pi(D)
       <=C sqrt(L_N[N^-2+eta_N]).                   (2)

The same pi(D) is used throughout; no unbounded preparation-TV transfer or
false uniform scaled-gradient bound occurs.

## 2. Auxiliary characteristics versus the paid SQRT ones

Write beta_Q^D,G_Q^D for the actual reflected SQRT characteristics from SQ.
On an internal edge P_app=Q_N a_N(p)a_N(p')/d_N. Uniform boundedness
and adjacent O(1/N) variation in AD-001 imply pointwise

    beta_app=a_N(p)^2 beta_Q^D+O(N^-1/4),
    G_app=a_N(p)^2 G_Q^D+O(N^-1).                    (3)

The errors follow by summing actual p_a R_ab^N and bounded sign jumps,
whose sums are uniformly bounded before critical acceleration. No parity
lumpability or fast-coordinate independence is used.

MT moments and |a_N-1|<=C(||p-u||^2+N^-1) give, with any fixed polynomial
Z weight, a_N^2-1=O(N^-1/2) in the required stationary Lp norms. The
already paid SQ remainder and boundary correction then give

    E_piD[(1+|Z|)|beta_app-f(Z)|]
       <=C[N^-1/4+N eta_N],
    E_piD|G_app-2q|<=C[N^-1/4+eta_N].               (4)

In particular the approximation was not engineered from f: its defining
equation came first from the full count scaling, and (3),(4) derive the
resulting critical behaviour from the actual rates/moments.

## 3. Final quantitative defect bounds (not a process error rate)

Combining (1)--(4),

    D_N:=E_piD[(1+|Z|)|beta_N-f(Z)|]
       <=C[N^-1/4 sqrt(L_N)+N^3/4 sqrt(L_N eta_N)+N eta_N] ->0,
    C_N:=E_piD|G_N-2q|
       <=C[N^-1/4+N^-1 sqrt(L_N)+sqrt(L_N eta_N)+eta_N] ->0. (5)

The constants are existence bounds under inherited MT premises, not
calibrated numerical finite-N error bars. These are the EXACT two payments
left open by IS. They concern the original unaugmented global source,
not a large-spectator limit or locally reselected entropy law.

The actual discrete predictable variance is G_N-h_N beta_N^2. Since the
first payment gives E|beta_N|=O(1), the exact bounded-jump estimate pays

    h_N E beta_N^2<=2N^-3/4 E|beta_N| ->0.            (6)

Thus (5) proves the centered bracket, not merely a raw quadratic mean.
The original actual-source critical path theorem follows in the next file.
