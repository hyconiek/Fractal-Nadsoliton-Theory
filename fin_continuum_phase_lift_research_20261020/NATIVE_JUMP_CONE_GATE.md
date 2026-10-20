# CPL-004: an exact native moment-cone gate that the continuum lift cannot skip

Date: 2026-10-20. The first identities are actual finite-count kinematics
for ANY single-copy generator, not assumptions about CPL coefficients.

## 1. Universal soft/fast trace identity

Let Z=N^(1/4)m and Y=sqrt(N)(p-p_m) in the Euclidean fast space. For a
count move a->b write delta=e_b-e_a and d_sigma=sigma_b-sigma_a. Then

    Delta Z=d_sigma N^-3/4,
    Delta Y=[delta-d_sigma sigma/12]/sqrt(N),
    ||Delta Y||^2=(2-d_sigma^2/12)/N.                (1)

Opposite parity gives 5/(3N); same parity gives 2/N. If q_cross,q_same
are the two TOTAL off-diagonal rate sums at the state, and time is the
KR critical clock sqrt(N), then EXACTLY

    G_N=Gamma_(sqrt(N)L)(Z,Z)=4q_cross/N,
    tr Gamma_(sqrt(N)L)(Y,Y)
       =(5/12)sqrt(N)G_N+(2/sqrt(N))q_same.          (2)

Departure rejection or zero diagnostic rates do not change the allowed
increments, so (2) holds at every actual state. A CN cap implies

    (5/12)sqrt(N)G_N<=tr Gamma_fast<=2C sqrt(N).     (3)

For the stationary KR class E G_N->2 is already paid. Thus the expected
fast bracket per unit critical time is at least (5/6+o(1))sqrt(N).
More strongly, for A_N=int_0^S G_N ds in[0,4CS], E A_N->2S. Eventually
P(A_N>=S)>=1/(8C). The integrated fast predictable bracket therefore
exceeds (5/12)S sqrt(N) with that positive probability and is NOT tight.

This rules out convergence of these full native joint martingale
CHARACTERISTICS to a fixed finite-bracket diffusion. It does NOT alone
rule out tight state paths with cancelling drifts. A two-state chain
{0,1/n}, flip rate n^3, has sup|X|<=1/n while its compensated martingale
bracket is nS. This is a precise logical countercontrol, not FIN.

## 2. Why the fixed-clock continuum process is not a native construction

The CPL process at kappa=1 has bounded local quadratic coefficients in
the original Y coordinates. Its full characteristics therefore cannot
be silently imported as the joint native critical characteristics in (2).
Neither a continuous phase construction nor a generic Markov approximation
theorem supplies actual count rates, all vertex drift equations, true E_N
conditions, fixed support or one common CN cap.

## 3. A scale-compatible comparison still is not a jump-cone solution

The comparison family itself can take kappa_N=sqrt(N) without changing
ANY law of Z or the nonaveraging identity. This only rescales its already
specified fast operator, not a physical clock or native count chain.
Its quadratic coefficients in whitened coordinates obey explicit bounds.
With p in[1/2,5/8], (p')^2<=1/32, and r=r_P(z,h),

    2a r_h^2<=40/9 in A, <=2 in B,
    2r^2/(1+r^2)<=2, r_z^2<=2/9.

These follow from -log(1-x)>=x. The pair's fast trace is at most 58/9;
the other eight OU directions contribute 16. The soft lift adds at most
(8/3)(2/9)=16/27. Unwhitening by fixed H_fast gives an upper bound C kappa_N.

Since H_rest<=12I, its eight OU directions alone contribute at least
(4/3)kappa_N in Euclidean Y trace. At kappa_N=sqrt(N), this exceeds the
necessary lower trace (5/12)sqrt(N)g_P by at least (2/9)sqrt(N), since
g_P<=8/3. Thus the comparison CAN meet the growth order and necessary
trace inequality; the fixed-clock mismatch is not an impossibility of
every scaled realization.

But TRACE compatibility is not membership in the full nonnegative
single-copy moment cone. The [cosine-variance witness](NATIVE_COSINE_VARIANCE_OBSTRUCTION.md)
actually REJECTS this frozen comparison's direct native characteristic transfer
on a set of positive stationary mass; it is not merely an unexecuted gate.
The actual matrices, drift at EVERY finite type,
Gibbs conductances, permitted phase/boundary edges and global capacity
must still be paid for any DIFFERENT prospective construction. No such cone/flow solution, native approximation or
finite-N IBA counterfamily is proved here. The original FP/RB/PC gates
remain open, rather than renamed or bypassed by this continuum object.
