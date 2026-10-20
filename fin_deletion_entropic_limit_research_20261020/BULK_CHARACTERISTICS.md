# DEL-002: actual bulk coefficients of the opposite-parity reference

Date: 2026-10-20. Statements are uniform on fixed compact subsets strictly
inside either limiting phase. No expansion at the deleted interface is used.

Let q=12, sigma_j=(-1)^j, E the fast space perpendicular to 1 and sigma,
and H=E_fast^T(qI-g_c A7)E_fast>0. The fixed-scaled compact and all its
single-copy neighbours lie in the true D_N\E_N eventually by MT. If it also
stays strictly inside its phase, all opposite-parity rates are retained.

## 1. Mean-field drift with the actual per-copy attempt normalization

Put s_j(p)=softmax_j(g_c A7 p), m=sigma.p, q_s=sigma.s. Dropping only
same-parity moves gives the exact mean-field formula

    B_j^cross(p)=[(1-sigma_j m)/2]s_j(p)
                         -p_j[(1-sigma_j q_s)/2].           (1)

The native LOO error is uniformly O(1/N), including its local derivatives.
On p_m=(1+m sigma)/q, (1) is sigma_j[tanh(m)-m]/q; its fast projection
is zero. For a fast perturbation w at the uniform point, both group-total
derivatives vanish and Ds(w)=g_c A7 w/q. Consequently

    DB_fast^cross(u)w=-(qI-g_c A7)w/(2q)=-H w/24.           (2)

On fixed (z,Y) compacts, m=N^-1/4 z and w=N^-1/2 Y. Smoothness, (2) and
the finite LOO error therefore imply L-dagger Y->-H Y/24.

## 2. Covariance factors and the zero mixed term

At uniformity in the leading mean-field coefficient, each ordered
opposite-parity count rate is N/144. Actual finite LOO rates have an O(1)
correction to that count rate; the already stated scaled errors charge it.
The complete bipartite K6,6 Laplacian acts by 6 on E, by 12 on sigma,
and by zero on 1. The projected jump frame thus gives ordinary fast
covariance I_E/12. Its mean on each oriented parity class is zero.

For A_N=sqrt(N)L-dagger the required compact limits are

    A_N Z -> b(z),                  Gamma_(A_N)(Z,Z)->2,
    A_N Y/sqrt(N) ->-H Y/24,        Gamma_(A_N)(Y,Y)/sqrt(N)->I/12,
    Gamma_(A_N)(Z,Y)->0.                                  (3)

The first line uses the actual UF cancellation, since removing same-parity
moves does not change Lm in the bulk. For the last line the prefactor of
the mean-field frame sum is N^1/4; that sum is EXACTLY zero on every parity
line, because the up/down projected frame means vanish separately. Its
off-line error is O(||w||+N^-1), so the mixed bracket is O(N^-1/4) locally.
There is no hidden O(N^1/4 m) mixed term.

These coefficients retain the true factor 1/2 in fast mobility. No clock
was refitted to erase the deletion of same-parity moves. In particular
the fast Gaussian generator is

    L_fast=(Delta_Y-(H Y).grad_Y)/24.

It is symmetric in gamma_H. This is only a local coefficient computation;
we do not infer full fast-path convergence at critical time, where the
fast predictable brackets diverge as shown by the prior native cone gate.
