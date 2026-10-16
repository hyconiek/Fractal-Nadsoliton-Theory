# SC-001 — bounded analytic attempt on the actual signed interface residual

**Disposition: OPEN_SIGNED_INTERFACE_LEMMA.** New proof: nine legal phase-safe parity-direction neighbours at every forbidden-edge tail, and an exact *cut-local* shielding bound. This does not settle (SC) for all genuine orbit-union cuts.

## 1. Exact frozen data, no invented graph

Let q=12, N in 12N, A=A7, g_c=12/Lambda6 exactly, D_N the genuine MT connected sublevel with rejection-reflected original LOO heat-bath, E_N its actual departure layer. Freeze the same functions and inclusive tie rule:

- sigma_j=(-1)^j; m(n)=N^{-1} sigma.n; z=N^{1/4}m;
- C(n)=sum_j n_j cos(pi j/2), S(n)=sum_j n_j sin(pi j/2), both integers;
- R=2(C^2+S^2)/N; gamma=1-Lambda3/Lambda6>0;
- p(z)=1/2+e^{-z^2}/8; T(z)=-(2/gamma) log(1-p(z)); A={R<=T}, B={R>T};
- U is the upward-oriented **opposite-parity** native graph inside D_N; H contains within-phase edges plus phase-crossing edges incident on E_N; F=U\H.

Write W_n=(N!/prod n_j!) exp[g_c n^T A n/(2N)], u_e=W_x q^0(x,y), s=A_U u, and for a G_+-orbit-union subset V:

`B(V)=sum_{x in V} W_x #{ H-upward edges x->y outside V }`, `Hout(V)=sum_{H-up V->Vc} u_e`, `Fout(V)=sum_{F-up V->Vc}u_e`, `Uin(V)=sum_{U-up Vc->V}u_e`.

The exact signed residual is `D(V)=Fout(V)-Uin(V)`; summing the native divergence yields `s(V)=Hout(V)+D(V)`. The required **new** all-N conclusion is `D(V)<=L N B(V)` on the unresolved orbit-union class R_N with one L independent of N. It remains OPEN. No rate modifications or new microscopic laws are used.

## 2. A genuine new finite q12 combinatorial fact: nine phase-safe moves

Put `Tmin=(2/gamma)log 2>0`, `Tmax=(2/gamma)log(8/3)`, `L_T=4/(3 gamma sqrt(2e))` and the previously paid shell width

`eta_N=4 sqrt(Tmax) N^(-1/2)+4/N+2 L_T N^(-3/4)`.

The inherited exact shell theorem implies every endpoint of every F edge lies in `|R-T(z)|<=eta_N`. Fix N sufficiently large that simultaneously

`eta_N<Tmin/2`, `sqrt(N*Tmin/4) >=1+(L_T/2)*N^(1/4)`.

Such an N0 exists because Tmin>0, eta_N=o(1), and the left side is order sqrt(N) while the growing term on the right is order N^(1/4). No claim about a practical N0 is made.

**Theorem (local nine-way escape).** If x is any endpoint of a forbidden F edge, then for EACH choice of parity direction (up and down) there exist **at least nine distinct legal neighbours in H** with that parity direction and in x's own phase. All of these nine move in the *radially phase-safe diagonal direction* in the integer (C,S) plane. Their native conductances obey `u_xy >= a0 N W_x`, `a0=(103/15000)e^{-24}`, using the already-paid MT native ellipticity input.

**Proof.** Every forbidden endpoint lies outside E_N, by the definition of H, so *every* one-copy proposal with a different target label is retained within D_N. The MT occupation floor ensures n_a>0 for all twelve possible departure labels. For upward moves, the departure odd label has Fourier vector `(0,+/-1)` and arrival even label has `(+/-1,0)`; hence every pair of signs `(Delta C,Delta S)=(+/-1,+/-1)` occurs in exactly 3x3=9 **different** odd-to-even label replacements. For downward moves, departure even and arrival odd again realize exactly nine replacements for each of the four sign combinations. These two directions change z by `+/- 2 N^(-3/4)` respectively.

Since `R>=Tmin-eta_N>Tmin/2`, we have `|C|+|S| >= sqrt(N*Tmin/4) >=1+(L_T/2)N^(1/4)`.

Choose a sign at zero arbitrarily. For the *inward* Fourier diagonal `(Delta C,Delta S)=(-sign C,-sign S)` the exact identity

`Delta R=4(C Delta C+S Delta S+1)/N =4(1-|C|-|S|)/N <= -2 L_T N^(-3/4)`.

For the *outward* diagonal `(+sign C,+sign S)`, `Delta R=4(1+|C|+|S|)/N>=2 L_T N^(-3/4)`. On either parity direction `|Delta T|<=2L_T N^(-3/4)` by the already paid Lipschitz bound. If x belongs to A, all nine inward-diagonal moves remain in A (including possible equality). If x belongs to B, all nine outward-diagonal moves remain strictly in B. Thus these are H edges. Native ellipticity on the genuine local domain yields the stated edgewise conductance bound. QED.

The result also gives the **pointwise** inequality at every upward forbidden tail `x`:

`sum_{F-up edges from x} u_e <= N W_x <= (9 a0)^(-1) sum_{H-up edges from x}u_e`.

This is genuine local microscopic compensation. It is **NOT** an inequality between outward H-cut edges and forbidden F-cut edges: the nine safe H heads may lie INSIDE an arbitrary cut V.

## 3. One actual cut subclass paid, and the remaining obstruction localized to shielded tails

For an orbit-union cut V put `T_F(V)={x in V: there exists a forbidden upward edge x->y with y outside V}`. Partition these tails into:

- `E_F(V)={x in T_F: there exists an upward H edge x->y with y outside V}` (exposed tails);
- `Z_F(V)=T_F\E_F` (shielded tails).

The native exit cap `sum_y q0(x,y)<=N` implies

`Fout(V)<= N sum_{x in E_F}W_x +N sum_{x in Z_F}W_x`.

Every exposed tail contributes at least one **distinct tail-weight unit** to `B(V)` (distinct tails count separately), so `sum_{x in E_F}W_x<=B(V)`. Therefore the following is a fully **proved finite q12 bound**:

`D(V) <= N B(V) +N W(Z_F(V)) - Uin(V)`, where `W(Z)=sum_{x in Z}W_x`.

Consequences: (i) for every V with `Z_F(V)=empty`, `(SC)` holds with `L=1`; (ii) if `(SC)` fails with some L>=1, then necessarily `N W(Z_F)-Uin >(L-1)N B`; (iii) each shielded tail has all its nine phase-safe H successors inside V, by the nine-way theorem. This is a concrete necessary geometric structure for a countercut on the *unchanged* q12 graph.

**Exact unpaid global signed absorption lemma:** for all genuine unresolved orbit-union cuts V and all sufficiently large N, bound

`N W(Z_F(V)) -Uin(V) <= L' N B(V)`

with one finite L', or supply a genuine diverging cut. Proving this would pay `(SC)` with L=1+L'. This is a strictly stronger sufficient strategy than (SC) itself and is not adopted as an axiom. No evidence from the finite controls establishes it. A pointwise escape property cannot imply it because the nine safe successors may all belong to V.

## 4. Why the proof stops here

Nothing in the shell mass bound, even `pi_D(Shell)->0`, controls the ratio `W(Z_F(V))/B(V)` when V concentrates near the shell. Nor do the exact pointwise native rate bounds force Uin to absorb shielded-tail source. To turn the nine moves into an outward cut charge one must prove a stopping/expansion property of the **real** H paths plus exact Gibbs ratio/capacity control on that path. Such a property was not proved, and it must not be postulated. The present proof is nontrivial local geometry and a restricted sufficient cut theorem, not a universal cut theorem.
