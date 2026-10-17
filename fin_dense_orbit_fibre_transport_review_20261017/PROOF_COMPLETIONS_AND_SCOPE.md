# DH dense-complement proof audit and the atom-mass payment

Date: 2026-10-17. Analytical completion under actual admitted MT/PC inputs;
not mechanized proof or peer review. Incoming source bytes are preserved.

## 1. The exact complement argument is valid

For A removed from the left fibre, set T=X_l\A and
R=X_l+1\N_H(T). Equal total currents give

    I'(N_H(T))-O'(T)=O'(A)-I'(R).

Two distinct count predecessors of one right vertex have L1 distance <=4.
Pairwise separation >4 in A therefore removes at most one predecessor of
any right state. Every nondeparture right state has at least nine legal
downward H predecessors, so R is contained in the true right E layer.
No within-phase replacement of full H or omission of left E rows is used.

The bulk row lower .173NW follows from .18NW incoming minus .007NW source.
For an E column, first bound the **uncorrected** I coefficient everywhere
on its phase fibre by .35NW (it is a fibrewise constant times W). Then
`I'_E=max(I,-s)<=NW`, since |s|<=NW. Do not infer that corrected E columns
satisfy the bulk .35 bound. The source uses the correct larger E bound.
Thus .173W(A)>=W(E_right) is a valid finite certificate once those current
bounds are paid at that N. They are eventual compact-band bounds, not a
usable numerical N0.

## 2. Explicit orbit and exact integer geometry

For the supplied ramp t, each parity sum is zero and its twelve G_+ images
are distinct, with exact minimum L1 distance36. Evenly distributing a
parity total gives entries floor(K/6) or ceil(K/6), so any two permuted
rounded bases differ in L1 by at most12. Consequently the count orbit's
minimum distance is at least36k_N-12. Components are positive for large N.

C and S equal -6k_N+O(1); R_N=O(N^-1/3)->0. The phase-A threshold is
bounded below by a fixed positive number. Thus all orbit members lie
strictly inside phase A eventually; no floating transcendental tie decision
is needed. With exact k_N=floor(N^(1/3)), the fast displacement is
O(N^-2/3), hence scaled Y=O(N^-1/6), while z stays in[-1,1].

The actual MT central connectivity includes these states and all their
one-copy neighbors, with F=O(1/N)<h/2 eventually. This puts A outside E,
not merely inside an invented cylinder.

## 3. Pay the polynomial atom mass directly

Uniform positive-interior Stirling bounds give, for the orbit representative,

    W(n*) >= c 12^N N^(-11/2),

because NF(n*/N) stays bounded uniformly on the compact soft band and
Y tends to zero. At fixed parity, the ten-dimensional fast lattice cell
volume is6N^-5. The inherited coercivity and bounded cell oscillations
bound its Gaussian lattice sum by C N^5, uniformly in lattice translation
and bounded z. Uniform Stirling then gives

    W(X_l) <= C 12^N N^(-1/2).

Therefore W(A)/W(X_l)>=c N^-5. This requires no numerical prefactor or
unpaid thin-shell limit. Uniform compact slice lower bounds also give
W(X_l+1)/W(X_l) bounded. The actual MT E-mass exponential estimate,
divided by polynomial slice mass, proves W(E_right)/W(X_l)<=poly(N)e^-bN.
The finite certificate in section1 follows uniformly for all sufficiently
large admissible N. The source's inherited atom estimate is thus paid by
the real current proof inputs, not by its finite orbit checks.

## 4. A computational qualification, not a theorem refutation

verify_dh.py chooses k through `floor(N**(1/3))`, and some fixture fibres
through a floating three-quarter power. Later orbit arithmetic is integer/
Fraction, but this is not an all-integer root API. For example, floating
evaluation at N=1200^3 gives1199 instead of1200 in this environment.
The mathematical definition and asymptotic proof are unaffected. The new
intake checker uses exact integer roots, checks perfect powers and rejects
an off-by-one control without editing the archived producer.

## 5. What remains genuinely open

Separated holes are not arbitrary clustered complements. Entire predecessor
clusters can be removed, producing right holes away from E. Small global
E mass alone does not pay Hall when a removed complement has even smaller
row mass. No full Hall, actual adverse family, band coupling, cap5N or
full-domain/IBA consequence follows. Do not default to more ramp patterns
or merely enlarge a separation constant.
