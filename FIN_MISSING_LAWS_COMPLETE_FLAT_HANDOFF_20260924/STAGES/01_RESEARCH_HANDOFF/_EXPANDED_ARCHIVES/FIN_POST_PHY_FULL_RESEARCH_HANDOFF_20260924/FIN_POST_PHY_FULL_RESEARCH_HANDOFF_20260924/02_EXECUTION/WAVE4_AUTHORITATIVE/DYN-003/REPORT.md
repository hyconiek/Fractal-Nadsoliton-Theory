# DYN-003 — Derived nonlinear phase potential and drift obstruction

Status: **PASS_CONTROLLED_PERIODIC_PHASE_REDUCTION_AND_EXACT_NO_DRIFT_FOR_UNBIASED_DISSIPATIVE_FRONTS**.

## Frozen inputs

The repository already proves all-orders relative-phase locking: positive-amplitude maxima satisfy `phi_k=k alpha`, and finite `q` reduces the common shift to a discrete `Z_q` orbit. Wave 3 strengthens this with an all-orders pinning bound: on bounded amplitudes the complete Umklapp barrier is exponentially small in `q`, while every fixed finite `q` with positive amplitudes remains pinned.

## Effective phase object

After transverse modes are eliminated inside the certified locking neighbourhood, the common phase has a `2pi/q`-periodic effective potential

`U_q(alpha)=u_0 + sum_(m>=1) [a_m cos(m q alpha)+b_m sin(m q alpha)]`,

with symmetry-related minima of equal energy. For the reflection-symmetric strict carrier the sine coefficients vanish after an origin choice. The PHA-002 analytic-strip estimate controls the total pinning amplitude rather than only the first Umklapp monomial:

`Delta_q <= 4 M_rho exp(-rho q)/(1-exp(-rho q))`.

A declared local phase action therefore has the form

`E[alpha]=integral [K/2 (partial_x alpha)^2 + U_q(alpha)] dx + controlled hard-mode corrections`,

where the corrections must be retained if they are not negligible. Symmetry-preserving variational corrections cannot create an energy difference between adjacent `Z_q`-related vacua.

## Exact no-drift identity

For overdamped gradient dynamics

`gamma alpha_t = K alpha_xx - U_q'(alpha)`

and a travelling front `alpha(x,t)=a(x-c t)` joining two minima, multiplication by `a'` and integration gives

`gamma c integral (a')^2 dz = U_q(alpha_+) - U_q(alpha_-)`.

Adjacent FIN translation vacua have equal energy, so if `gamma>0`, every finite-energy travelling front must have **c=0**. The same identity holds for damped inertial dynamics

`m alpha_tt + gamma alpha_t = K alpha_xx - U_q'(alpha)`.

Inertia changes the profile equation but not the integrated balance: persistent dissipative drift still requires a nonzero well-energy difference or a nonconservative drive.

For purely conservative dynamics (`gamma=0`) equal-energy kinks may move if a suitable inertial model is separately supplied; their velocity is then an initial-condition/dynamical parameter, not an autonomous FIN drift source.

## Joint q,N interpretation

Two limiting regimes are now sharply separated:

- fixed finite `q`: the common phase is pinned (`Delta_q>0`), so an unbiased gradient front is static and lattice/Umklapp pinning remains;
- `q` growing so that `N Delta_q ->0`: pinning can disappear, but then the localization barrier supporting a discrete domain wall also vanishes unless another sourced nonlinear term survives.

Thus the present phase reduction does **not** by itself provide an autonomously moving particle-like excitation. A bias, drive, separately sourced inertia with prepared momentum, or other non-gradient kinetic ingredient must be declared.

## Held-out controls

Reversing an explicit bias reverses the sign of the energy difference and therefore the front speed sign. The no-drift result is independent of a fitted `q`: it applies to every finite `q` for which the effective potential has exactly degenerate translation-related minima.

## Allowed conclusion

FIN supplies a controlled periodic phase-pinning structure and, under declared dissipative dynamics, an exact transport obstruction. It does not yet supply a pump, physical inertia, stable moving soliton, or matter ontology.
