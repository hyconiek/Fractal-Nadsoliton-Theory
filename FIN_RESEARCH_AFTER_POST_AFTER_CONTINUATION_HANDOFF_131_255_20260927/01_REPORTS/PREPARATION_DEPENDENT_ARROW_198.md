# PREPARATION-DEPENDENT-ARROW-198
## The reversible FIN heat bath has no stationary arrow; monotone relaxation appears only after nonequilibrium preparation

Date: 2026-09-26

Status:
exact information-theoretic theorem plus exact N=6 numerical replay.

## 1. Equilibrium path law

The exact leave-one-out Gibbs generator satisfies detailed balance:

    pi_i L_ij
      =
    pi_j L_ji.

For N=6 at

    g=5.145228719489142

the maximum numerical detailed-balance residual in the reconstructed exact
generator is

    1.69e-17.

Hence the stationary path law is time-reversal invariant.

The stationary entropy production rate is zero.

So the microscopic law itself selects no thermodynamic future direction.

## 2. Nonequilibrium preparation

Let p_t be any distribution evolving under the Markov semigroup P_t, with
stationary pi.

Relative entropy obeys the data-processing inequality:

    D(p_t P_s || pi P_s)
      <=
    D(p_t || pi).

Since

    pi P_s=pi,

    boxed:
    D(p_(t+s)||pi)
      <=
    D(p_t||pi).

Thus nonequilibrium free-information distance to equilibrium is monotone.

This creates an operational relaxation arrow once an initial preparation is
specified.

## 3. Exact N=6 replay

Starting from the extreme count state

    all six copies on label 0,

the exact KL divergence decreases:


    t=0:
      D(p_t||pi)=2.886663851349

    t=0.02:
      D(p_t||pi)=2.859127424010

    t=0.05:
      D(p_t||pi)=2.830792179796

    t=0.1:
      D(p_t||pi)=2.794590433677

    t=0.2:
      D(p_t||pi)=2.741460111208

    t=0.5:
      D(p_t||pi)=2.644993213390

    t=1:
      D(p_t||pi)=2.557381050070

    t=2:
      D(p_t||pi)=2.443705709003

    t=4:
      D(p_t||pi)=2.228997141205

    t=8:
      D(p_t||pi)=1.824763829013

    t=16:
      D(p_t||pi)=1.253014024469


No reversal-breaking term was added to the generator.

The asymmetry resides in the boundary condition:

    special prepared p_0
      ->
    ordinary equilibrium future.

## 4. Global reversible dilation

Report 188 embeds the same observed process in an invertible two-sided history
shift.

There the complete microscopic-history law remains reversible.

The observed relaxation arrow therefore comes from:
- conditioning on a special preparation;
- discarding the full hidden past/future history.

This is consistent with global information continuity.

## 5. Conclusion

FIN now distinguishes cleanly:

    temporal orientation / history coordinate
      from
    thermodynamic arrow.

The first follows from choosing an orientation of the transformation chain.

The second requires a nonequilibrium preparation or equivalent low-entropy
boundary condition.

Neither fixes seconds.
