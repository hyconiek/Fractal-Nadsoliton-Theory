# FP-001 — exact original MT domain, fixed phase predicate and actual drift source

This is a source-contract audit against accessible older **source** MT-001/003, BA-001 and KR-002 files. The later BA review `COUNT_PARTITION_FEASIBILITY_SEED.md`, its review/verification, and full current repository AGENTS are unavailable; they are not claimed to have been read.

For `p=n/N`, `sigma_j=(-1)^j`, `u_j=1/12`, `m=sigma·p`, `p_m=u+m*sigma/12`, `w=p-p_m`, the actual MT analysis cylinder uses `|m|<.01`, `||w||_2<.0001`, all coordinates positive, and `F_gc(p)<=1e-10`. Its count domain `D_N` is the nearest-neighbour component containing the uniform count vector. MT-001 proves that the full sublevel is connected for all sufficiently large multiples of 12, conditional on its admitted exact spectral gap.

`F_g(p)=sum_j p_j log(12 p_j)-(g/2)*p^T A7 p`; the Gibbs weight is `W_N(n)=N!/(prod_j n_j!)*exp[g_c/(2N)*n^T A7 n]` on `D_N`, where inessential constants `12^{-N}` and normalization are suppressed. Never substitute the global equilibrium, artificial enlarged domain or a rounded numerical gain.

`C(n)=sum_j n_j cos(pi*j/2)` and `S(n)=sum_j n_j sin(pi*j/2)` are **integers** (the original k3 Fourier sums). Write `Z=N^(1/4)m`, `R=2(C²+S²)/N`, `p(Z)=1/2+exp(-Z²)/8`, `gamma3=1-Lambda3/Lambda6` and `T=-(2/gamma3)log(1-p)`. Freeze `A={R<=T}`, including ties, and `B=D_N\A`. The exact equivalent predicate is `exp(-gamma3 R/2)>=1/2-exp(-Z²)/8`, so numerical evaluation must return INDETERMINATE if rigorous interval bounds overlap. This equivalence alone does not provide an exact algebraic decision procedure for all spectral/transcendental comparisons.

Let `a` be the departed label and `b!=a` its new label. The original reflected count rate **only for n,n-e_a+e_b both in D_N** is

`q0(n,n-e_a+e_b)=n_a * exp[(g_c/N)(A7(n-e_a))_b]/sum_j exp[(g_c/N)(A7(n-e_a))_j]`.

Leaving D_N is rejected, not renormalized. Define exact `d0(n)=sum_b q0(n,n-e_a+e_b)*(sigma_b-sigma_a)/N`, summing all admissible `a,b`, and **the full statewise demand** `b_n=W_N(n)d0(n)` (or equivalently `pi_n d0(n)`). No replacement by a cubic approximation or a conditional mean is allowed in the finite primal.

Cyclic rotation by one label sends `m->-m`, `(C,S)->` a right-angle rotation, leaving `R`, `T`, `A/B`, `D_N`, Gibbs weight invariant. Native drift is anti-equivariant. Therefore `sum_{n in A} W_N(n)d0(n)=sum_{n in B} W_N(n)d0(n)=0` **exactly**, a necessary phase-total balance. It is not sufficient for pointwise feasibility.

MT-003 inherited conditional fast Gaussian with `Var(C)=Var(S)≈N/(2 gamma3)` and independence from the quartic slow coordinate yields `R=>chi²_2/gamma3`; its CDF is `1-exp(-gamma3 r/2)`. Thus the *already submitted* static phase mass function is `P(A|Z=z)=p(z)`. This is an asymptotic equilibrium statement, not existence of closed-phase rates.
