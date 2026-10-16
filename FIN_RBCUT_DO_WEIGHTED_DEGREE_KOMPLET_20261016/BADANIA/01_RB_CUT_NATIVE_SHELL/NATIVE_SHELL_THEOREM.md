# RB-CUT-SHELL: exact Fourier step geometry and concentration of all forbidden edges

**Status:** NEW EXACT FINITE-N STRUCTURE, followed by a **conditional** stationary-asymptotic consequence using MT-003. **NOT** an exact lower-ideal sign proof or legal diagnostic generator.

## Frozen genuine FIN objects

Retain original q12/A7, exact `g_c=12/Lambda6`, `N` a multiple of 12, true MT-001 local count sublevel `D_N` and `E_N` as its genuine native departure layer. Set

- `sigma_j=(-1)^j`, `m(n)=sum_j sigma_j n_j/N`;
- `C(n)=sum_j n_j cos(pi j/2)`, `S(n)=sum_j n_j sin(pi j/2)`; both integers;
- `Z_N(n)=N^(1/4)m(n)`, `R_N(n)=2(C(n)^2+S(n)^2)/N`;
- `gamma_3=1-Lambda3/Lambda6>0`;
- `p(z)=1/2+(1/8)exp(-z^2)`, `T(z)=-2 log(1-p(z))/gamma_3`;
- `A_N={n in D_N:R_N(n)<=T(Z_N(n))}`, **inclusive ties**; `B_N=D_N\A_N`.

`U_N` comprises all native opposite-parity one-copy edges inside `D_N`. `H_N` retains all within-phase edges and any cross-phase edge incident on `E_N`; `F_N=U_N\H_N` are the forbidden interior cross-phase edges. These are not physical domain exits.

## Theorem 1: finite-N forbidden-edge annulus

For every opposite-parity replacement `a -> b`, one index has even and one has odd parity. The vectors `(cos(pi j/2),sin(pi j/2))` are coordinate axes with unit lengths. Thus

`Delta C=+/-1`, `Delta S=+/-1`, `Delta C^2+Delta S^2=2` (**exact integers**), and

`R_N(y)-R_N(x)=[4(C_x Delta C+S_x Delta S)+4]/N`, `|Z_N(y)-Z_N(x)|=2 N^(-3/4)`.

Put

`T_* = (2/gamma_3) log(8/3)`,
`L_* = 4 / (3 gamma_3 sqrt(2e))`,
`eta_N=4 sqrt(T_*) / sqrt(N) + 4/N + 2 L_* N^(-3/4)`,
`Shell_N={n in D_N:|R_N(n)-T(Z_N(n))|<=eta_N}`.

The bounds `0<T(z)<=T_*` and `|T'(z)|<=L_*` hold for all real z. (For the latter, `T'(z)=(-z exp(-z^2)/(2 gamma_3))/(1/2-exp(-z^2)/8)` and `sup |z|exp(-z^2)=1/sqrt(2e)`.)

Take any phase-crossing edge `{a,b}` with `a in A_N`, `b in B_N`. Since `R(a)<=T(Z(a))<=T_*`, the **exact** Fourier jump identity gives

`|R(b)-R(a)| <= 4 sqrt(T_*)/sqrt(N)+4/N`.

By the Lipschitz bound, `|T(Z(b))-T(Z(a))|<=2 L_* N^(-3/4)`. Combining `R(a)<=T(Z(a))` and `R(b)>T(Z(b))` then yields **both endpoint inequalities**:

`0 <= T(Z(a))-R(a) <= eta_N`, `0 < R(b)-T(Z(b)) <= eta_N`.

Therefore **EVERY forbidden edge is incident only on vertices of `Shell_N`**. In fact this holds even for allowed cross-phase edges near `E_N`; the distinction between their permissions is preserved. The theorem applies to every genuine count vector and every true MT domain edge; it uses no numerical spectral rounding and does not assume a diagnostic flow exists.

## Theorem 2: quantitative native-flow localization and cut deficit

Let `W_x=(N!/prod_j n_j!) exp[(g_c/(2N)) n^T A7 n]`, `Z^W_N=sum_{x in D_N} W_x`, `pi_x=W_x/Z^W_N`, and `c^0_xy=W_x q^0(x,y)=W_y q^0(y,x)` for native reflected LOO heat bath.

At each count state the total labelled attempt rate is `<=N`, including after reflection. Orient each opposite-parity edge upward in `m`. Since both endpoints of any `F_N` edge lie in the shell, charging each forbidden edge to its lower endpoint gives

`sum_{e in F_N} c^0_e <= N sum_{x in Shell_N} W_x = N Z^W_N pi_D(Shell_N)`.

Define `s_x=(N/2)W_x L^0_D m(x)=sum_{e upward out of x}c^0_e-sum_{e upward into x}c^0_e`. For every H-lower ideal `L`, RB-001's exact identity implies

`sum_{x in L}s_x >= - sum_{e in F_N}c^0_e >= - N Z^W_N pi_D(Shell_N)`.

Equivalently the **most negative normalized native ideal source** is at least `-pi_D(Shell_N)`. No sign is asserted at finite N.

For another useful consequence retain exactly native conductances on `H_N` and delete all `F_N` conductances (same-parity native edges may remain). This creates an actual finite reversible *comparison* generator `L^H_N` with invariant `pi_D`, support in `H_N`, and total rate `<=N`. It is **not** in the frozen KR/RB diagnostic class because its statewise parity drift need not equal the native drift. The residual is shell-supported and obeys

`|L^H_N m(x)-L^0_D m(x)| <= 2 1_Shell(x)`;
`E_pi |L^H_N m-L^0_D m| <= 2 pi_D(Shell_N)`;
`0 <= G^0_N(x)-G^H_N(x) <=4 1_Shell(x)` where `G=N Gamma(m,m)`;
`E_pi (G^0_N-G^H_N) <=4 pi_D(Shell_N)`;
`||s-A_H c^0_H||_1 <=2 sum_{e in F_N}c^0_e <=2N Z^W_N pi_D(Shell_N)`.

The `G` inequality is exact because opposite-parity jumps have squared increment `4/N^2`; same-parity jumps have zero increment. **An approximately drift-matched generator is not a legal exact-drift flow.**

## Theorem 3: stationary asymptotic shell mass (inherited MT premise)

Under the already admitted stationary local MT-003 critical/Gaussian law (which is a separate scientific predecessor), `(Z_N,R_N)` converges in distribution to `(Z,R)` with `Z` quartic and, conditional on Z, `R` having the continuous `chi^2_2/gamma_3` law. Thus the boundary `R=T(Z)` has probability zero. Since `eta_N->0` and the threshold is continuous,

`pi_D(Shell_N) -> 0`.

Consequently, uniformly over **ALL** (even very irregular) H-lower ideals at each N,

`inf_L [sum_{x in L}s_x/(N Z^W_N)] >= -o(1)`.

This is a new **normalized asymptotic defect bound**, NOT the necessary exact inequality `sum_L s>=0`. MT gives no useful local anti-concentration rate for this moving shell from its weak convergence alone. Even `o(N)` unnormalized-by-probability scale can be negative. Neither actual routing nor a uniform cap is proved.

## Irreducible remaining dependency

To solve RB-CUT, one must **either** establish the exact nonnegative sign for all actual `H_N` ideals (not only that the minimum tends to zero after normalization), **or** supply exact negative ideals with valid N quantifiers. If the sign is proved, uniform vertex capacity is still an independent gate. The theorem does not validate any AS/KB physical source or ToE proposition.
