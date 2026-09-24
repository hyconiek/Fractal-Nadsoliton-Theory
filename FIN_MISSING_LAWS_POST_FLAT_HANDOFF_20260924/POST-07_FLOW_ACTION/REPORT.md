# POST-07 — path-space action no-go

Status: **NO_GO_EQUILIBRIUM_DOES_NOT_SELECT_PATH_SPACE_ACTION**

Detailed balance fixes only rate ratios.  A general reversible family can be
written

`k_ij = a_ij sqrt(pi_j/pi_i)`, `a_ij=a_ji>0`.

The symmetric activities are kinetic data not contained in pi.

Using the repository's local q=12 test with one frozen neighbor:
- Barker/Metropolis conductance ratio ranges from
  `0.500000` to `0.621403`;
- continuous-time spectral gaps are
  Metropolis `0.799818`,
  Barker `0.467185`;
- stationary path relative-entropy rates are
  M||B `0.154293`,
  B||M `0.124650`.

The nonconstant conductance ratio proves these dynamics are not related by a
single time rescaling.  Equilibrium FIN therefore does not source a unique
path-space action.
