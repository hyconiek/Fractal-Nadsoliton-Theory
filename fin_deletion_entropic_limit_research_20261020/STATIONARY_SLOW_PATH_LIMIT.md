# DEL-004: the actual stationary soft/phase process limit of L-dagger

Date: 2026-10-20. Conditional analytical theorem for the existing comparison
generator, NOT for primary HB or a drift-matched KR family.

## 1. Statement

Start the finite count chain X_N exactly from pi_D,N, unpostselected. For
fixed S let

    Z_N(t)=N^(1/4)m(X_N(sqrt(N)t)), P_N(t)=phase(X_N(sqrt(N)t)).

Then on D([0,S],R times {A,B}), with the usual J1 topology,

    (Z_N,P_N) =>(Z,P),
    P(t)=P0 constant,
    dZ=[-Z^3/3+partial_z log p_(P0)(Z)]dt+sqrt(2)dB.        (1)

The initial joint law is q_P(z)dz, so Pr(P0=P)=mu_P and
Z0|P0=P has density q_P/mu_P. Future Brownian increments are independent
of this initial pair in the natural Markov realization. No full fast-path
limit or inherited full count-filtration Brownian identification is claimed.

## 2. Finite-dimensional distributions use the FULL finite semigroup

DEL-003(7) applies to bounded smooth phase/soft sections. Bounded
multiplication by another such section passes through the cell embedding:
the piecewise constant multiplier converges off a null phase boundary and
is uniformly bounded. Iterating the full count Markov semigroup therefore
gives every stationary multi-time moment of bounded slow sections in the
limit. The prelimit observed (Z_N,P_N) is NOT assumed Markov.

The weighted form in (6) has a conservative one-dimensional diffusion in
each component. Its drift has the cubic inward part plus a bounded smooth
correction, since p_P>=3/8 and p'_P is bounded. Local Lipschitz coefficients
and the z^2 Lyapunov bound give nonexplosion/uniqueness; usual stopped
continuous dependence gives a Feller version on each line. This identifies
the limiting stationary finite-dimensional law as (1), not just a formal
stationary density. The two-component slow Hilbert space is the range of
intrinsic phase/soft functions, not an appended hidden state.

## 3. Tightness avoids an unproved temporal boundary-drift replacement

The critical soft jump is at most 2N^-3/4 and its martingale bracket density
is G_N<=4. Let M_N be the forward compensated martingale for Z_N on [0,S],
and Mhat_N that of the stationary reversed finite count chain. Reversibility
gives the same generator and the same deterministic bracket bound. At
every deterministic nonjump t,

    Z_N(t)-Z_N(0)
       =[M_N(t)-(Mhat_N(S)-Mhat_N(S-t))]/2.                 (2)

Use simultaneous rational times and the actual left/right reverse versions;
the pathwise endpoint discrepancy is at most a fixed multiple of the
vanishing maximal jump. The sign is MINUS for state increments, unlike
the PLUS identity for additive Poisson functionals.

Both martingales are C-tight by the bounded bracket, optional isometry
and vanishing jumps. Reversal is continuous at their continuous limits;
their joint tightness and (2), together with the MT initial moments, give
C-tightness of Z_N. This also bounds maximum containment, without
estimating the large interface drift by its absolute value in time.

Stationary phase-crossing intensity gives

    Pr(any P_N change before S)<=2S N^(3/2)pi_N(E_N)->0.

Thus the pair is tight with a constant-phase limit. Its already identified
finite-dimensional laws uniquely give (1). No pointwise convergence of
beta_N, unknown full-chain gap or assumed iid observations entered the proof.

## 4. What the theorem leaves open

The comparison's instantaneous noise DOES satisfy NC, but its drift is not
the native statewise drift: (1)'s additional force agrees with the independent
weak boundary profile (3). Therefore this is NOT a counterexample to KR
under its actual exact-drift hypotheses and does not decide universal native
IBA. Full drift-matched phase flows/Hall/gluing are still missing.

The theorem is stationary, reflected-domain and fixed-horizon only. No
late-iid transfer, unreflected comparison theorem, numerical N0/error rate,
all-time mixing or physical source is exported. Strong temporal replacement
of the raw interface drift, joint fast trajectories and entropy convergence
are not asserted merely from this weak path convergence.
