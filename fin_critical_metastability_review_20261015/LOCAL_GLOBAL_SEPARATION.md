# Local equilibration and a one-sided global mixing lower bound

Date: 2026-10-15. Status: analytical corollary of the admitted MT and UF results.

MT admits polynomial local equilibration before exit up to exp(cN), for
each fixed `0<c<b=5*10^(-11)` and the stated N-subsequence/initial ensemble.
UF's exact exclusion certificate applies eventually to the same gains.

The MT cylinder satisfies

    ||p-u||^2=m^2/12+||w||^2
       < epsilon^2/12+delta^2 < (1/10)^2.

Thus D_N lies in the excluded ball, and full Gibbs mass of D_N is at most
`(N+1)^12 exp(-11N/20)`. Actual iid-start trapping gives

    (mu0_N P_exp(cN))(D_N)>=1-o(1).

Therefore

    TV(mu0_N P_exp(cN), pi_global,N)>=1-o(1).

TV of a mixture is no greater than the worst starting-state TV. Hence the
**worst-case count-chain mixing time** to any fixed threshold below1, for
example1/4, exceeds exp(cN) for all sufficiently large admissible N.

This is a valid global **lower bound only**, alongside a polynomial **local**
relaxation guarantee. It does not identify the true global communication
height, supply a matching upper bound, solve Gamma=B4, compute capacities or
give matching residence/escape-time asymptotics. It does not calibrate a clock.

The source's broad "no global mixing" nonconclusion should be read with this
explicitly derived one-sided corollary, not used to hide it or to claim a sharp
global result. Local/sector equilibrium remains distinct from full Gibbs.
