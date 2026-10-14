# UF-001 — Exact reference point, unbiased preparation and spectral contract

## Model and exact parameter

The mathematical parent has q=12 labels, the *exact-decimal* kernel parameters omega=0.18575, phi=0.16250, eta=1.8, and the mathematical pi, with A=XX^T from the earlier rank-seven Fourier construction. The normalized count p=n/N evolves under the **leave-one-out**, rate-one-per-labelled-copy generator L_N; its finite-N Gibbs measure and the diagnostic MET/BAR rules are as accepted in KI. No empirical-refresh finite-N substitution is allowed.

Let u=(1/12,...,1/12), sigma_j=(-1)^j, m=sigma^T p, p_m=u+m sigma/12, lambda_6 the simple largest Fourier eigenvalue of A, chi=lambda_6/12, **g_c=1/chi defined by the exact spectrum**, g_N=g_c+lambda N^(-1/2). Numeric floating g_c is *not* substituted into a sequence N->infinity. The specific report's nominal numerical proposal is in DIAGNOSTICS.json. The accepted earlier review bounds [5.12342755139859,5.12342755139864] refer to the same exact-decimal spectral construction but this package does not reproduce that original interval source file.

For q=12 circulant A, J=g_c Sigma_u A-I on the normalized tangent, Sigma_u=(1/12)I on T={x:sum x=0}. The parity eigenvector sigma has zero eigenvalue, and every other one of the ten tangent Fourier directions has eigenvalue g_c lambda_k/12-1<0, using the *accepted spectral boxes* and unique-top-mode proof, including four absent k=1,2 directions with eigenvalue -1. The worst fast rate from the accepted earlier box is at least 0.01860477465101. Let R_f=I-11^T/12-sigma sigma^T/12 be the orthogonal fast projector. This is an analysis projection in the occupation simplex, not a new physical component. The exactly invariant parity line is the fixed space of the label-shift-by-two subgroup.

## Unbiased finite-N start

Take N iid labels with P(label=j)=1/12 independently of every kinetic realization. If p_N(0) is their empirical occupation, w_N(0)=R_f(p_N(0)-u), then bounded independent centered summands imply

E ||w_N(0)||^2 <= C_2/N, E ||w_N(0)||^4 <= C_4/N^2.

Proof of fourth moment: in each coordinate, expand (sum_{k=1}^N (1{Y_k=j}-1/12))^4; terms with an index appearing once vanish, while all remaining terms have O(N^2) count and bounded magnitude. Sum over the fixed 12 coordinates and apply (sum x_j^2)^2 <=12 sum x_j^4. Similarly m_N(0)=N^-1 sum sigma(Y_k) has variance 1/N and E m_N^4=(3N^2-2N)/N^4, hence N^(1/4)m_N(0)->0 in L2, with no seeded sign.

A balanced deterministic count vector at u (N divisible by 12) is a distinct comparison initial law and is never substituted for these multinomial starts.

## Stops and precise order of limits

For fixed R,S<infinity and bounded lambda, define Z_N(s)=N^(1/4)m_N(N^(1/2)s), a soft first-hitting time tau_{N,R} for |Z_N|>=R, fast tube tau_{N,f} for ||w_N||>=N^(-1/5), and horizon S. Freeze Z at its first soft boundary hit and record the capped time and side; the fast auxiliary stop is removed only after its probability is shown to vanish. This is a **mathematical critical gain sequence**, not a laboratory ramp or a source of cosmological initialization.

## Deduplication and status

Known earlier: parity eigenmode and pitchfork, exact equilibrium landscape, cubic mean-field drift on p_m, and non-Markov finite-N parity readout. New target: path-level weak convergence **from iid uniform actual label preparations**, a ten-fast-direction energy bound over physical model time N^(1/2)S, and the comparison with full Gibbs equilibrium. Primary microscopic kinetics is supplied, not derived by the theorem. Existing finite-N nonlumpability is consistent with asymptotic closure.
