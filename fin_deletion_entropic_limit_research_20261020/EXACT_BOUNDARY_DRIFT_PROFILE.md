# DEL-001: the actual stationary interface drift profile

Date: 2026-10-20. This theorem does not require a process limit.

Write pi_N=pi_D,N, A_N=sqrt(N)L-dagger, beta_N=A_N Z=N^(3/4)L-dagger m,
G_N=Gamma_(A_N)(Z,Z), and beta0_N=N^(3/4)L0_D m. Let rho be the normalized
exp(-z^4/12) density, b=-z^3/3, q_P=rho p_P, P=A or B.

## 1. Inherited facts are used in their actual scopes

RBCUT already proves G_N->2 in stationary L1: deleting same-parity moves
does not affect the soft bracket; deleted opposite-parity edges lie in the
Fourier interface shell of vanishing mass, and 0<=G0-G_N<=4 1_shell.
The actual rate cap gives 0<=G_N<=4 and soft jumps at most 2N^-3/4.

Every remaining phase change has an E_N endpoint. Reversibility and the
rate cap give stationary phase-crossing intensity at most 2N pi_N(E_N).
It follows on the critical clock that

    |E_pi[1_P A_N F]|=|E_pi[F A_N 1_P]|
                   <=2||F||_infinity N^(3/2)pi_N(E_N)        (1)

for every bounded count function F. This tends to zero exponentially.
Phase invariance is not assumed exactly. The Fourier shell is not E_N.

## 2. Test by a bounded antiderivative, not a compact test of the wrong order

For h in C_c^2(R), choose f(z)=integral_(-infinity)^z h(u)du. This f is
bounded, f'=h and f'''=h''. On every finite count space it is a legitimate
generator test. Taylor expansion and the cap give

    A_N f(Z)=beta_N h(Z)+(G_N/2)h'(Z)+R_N,
    ||R_N||_infinity<=C_h N^-3/4.

Use (1), G_N->2 in L1 and the actual MT joint phase/static limit. Then

    E_pi[1_P h(Z_N) beta_N]
              ->-integral q_P(z)h'(z)dz
               =integral q'_P(z)h(z)dz.                    (2)

This identifies the drift as a DISTRIBUTION against all compact smooth h.
No total-variation compactness of signed drift measures or pointwise/L1
convergence of beta_N is asserted. The bounded antiderivative permits
nonzero integral h, so no undetermined constant current is left over.

The native weighted drift theorem supplies beta0_N->b(Z_N) in stationary
weighted L1. Subtracting it from (2) yields the explicit new profile

    E_pi[1_P h(Z_N)(beta_N-beta0_N)]
                  ->integral h(z)rho(z)p'_P(z)dz.           (3)

Thus the effective conditional candidate drift is

    b_P(z)=q'_P/q_P=b(z)+partial_z log p_P(z).              (4)

This is a boundary/conditional-volume correction in the given comparison,
not an independently sourced physical force.

## 3. Checks against the older obstruction and the global cancellation

Summing (3) over both phases gives zero for every h, since p'_A+p'_B=0.
So globally averaged compact soft drift tests can miss the entire effect.
Phasewise signed tests do not. Dual approximation of sign(p'_P) gives

    liminf E_pi|beta_N-beta0_N|
                   >=(1/2)E_rho[|Z|exp(-Z^2)]>0.

The weighted version gives the ALREADY KNOWN RBCUT lower bound M/2,
M=E_rho[Z^2exp(-Z^2)], with signed phase contributions -M/4 and +M/4.
That old bound is a consistency control, not a new theorem counted twice.
None of these static identities alone pays temporal drift replacement;
the separate form/semigroup proof below supplies the process law.
