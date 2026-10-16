# A bounded proof tool for the next IBA gate

Date:2026-10-16. Self-contained finite reversible-chain argument. Standard
Poisson/forward-backward martingale tools are reused, not claimed as new
general probability theory. The unpaid FIN application is identified below.

## 1. Separate invariant components before taking an inverse

Let L_N be a finite generator reversible with positive pi_N. The KR model
class has no explicit irreducibility/lower-rate assumption. The connectedness
of the MT domain graph alone does not imply every diagnostic generator uses
all its edges. Let P_N be orthogonal projection in L2(pi_N) onto ker L_N:
it is conditional expectation on the actual communicating component.

Put r_N=G_N-E_pi G_N, h_N=P_N r_N and r_N^perp=r_N-h_N. The invariant part
stays constant along each path. Jensen at the final time gives

    E sup_(t<=S)|integral_0^t r_N(X_(sqrt(N)s))ds|
       >=S E|h_N|.

Thus vanishing component means are necessary for IBA. Global centering does
not suffice. For irreducible L_N, h_N=0 because r_N is globally centered.

## 2. Exact Poisson/Dirichlet quantity

Solve -L_N u_N=r_N^perp with zero mean on each component, and define

    R_N=<r_N^perp,u_N>_pi
       =D_N(u_N,u_N),
    D_N(f,f)=1/2 sum_(x,y) pi_x q_xy(f_y-f_x)^2.

The inverse is on the orthogonal complement of the kernel, not an assumed
full-space inverse. With finite positive rates on each component this
quantity is finite. It also equals

    sup_f {2<r_N^perp,f>_pi-D_N(f,f)}.

Complete the square with u_N to prove the identity; component constants
do not change it. Replacing r_N^perp by r_N makes that supremum infinite
when h_N!=0. No uniform spectral gap has been assumed or proved.

## 3. A maximal bound with the actual critical time factor

In scaled time A_N=sqrt(N)L_N, let v_N=u_N/sqrt(N). Then -A_N v_N=r_N^perp,
and E Gamma_A(v_N)=2R_N/sqrt(N). For the stationary process on [0,S] form
the forward martingale M and the martingale Mhat of the reversed process.
Reversibility makes its generator the same. At every deterministic t,

    integral_0^t r_N^perp(X_(sqrt(N)s))ds
       =[M_t+Mhat_S-Mhat_(S-t)]/2.

There is almost surely no jump at a specified deterministic time. Taking
the identity on rational times, and using continuity of the integral,
avoids ambiguity about reversed left/right jump conventions in the supremum.
Both martingales have second moment2S R_N/sqrt(N) at the final time.
Doob's inequality and |Mhat_S-Mhat_(S-t)|<=2sup|Mhat| give

    E sup_(t<=S)|integral_0^t r_N^perp ds|^2
       <=20S R_N/sqrt(N).

Including the component term gives the safe sufficient bound

    E sup_(t<=S)|integral_0^t r_N ds|^2
       <=2S^2 ||h_N||_pi^2+40S R_N/sqrt(N).

For irreducible or component-centered data the direct20 bound applies.
Constants are not optimized. Therefore

    ||h_N||_pi ->0 and R_N=o(sqrt(N))                 (PC)

imply IBA and hence KR's conditional process theorem. The deterministic
KR bound on G also bounds the additive functional on each fixed horizon;
the sufficient L2 estimate thus entails its required L1 estimate.

## 4. What this does not establish

PC is another sufficient route, not a necessary criterion or an automatic
consequence of NC. It can pay time averaging when pointwise NC fails;
conversely very small but extremely slowly moving r_N may satisfy NC while
failing this sufficient bound. Do not identify failure of PC with failure
of IBA or different critical dynamics.

The exact ten-state Poisson fixture in verification.json checks the equation,
energy and variational identities; it does not bound R_N for any q12
critical sequence. An abstract disconnected-state control illustrates the
component issue, not an admissible FIN counterexample.

The real open work is to estimate h_N and R_N from **actual KR count
conductances, exact statewise parity drift and pi_D,N**, or build a fully
admissible family with a genuinely nonvanishing bracket functional.
A native-HB gap, static conditional noise identity, large R_N alone or a
slow synthetic hidden variable does not pay that work. No physical source
or additional noise law is introduced by this proof tool.
