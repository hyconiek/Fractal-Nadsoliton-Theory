# MEMORY-RENORMALIZATION-OF-12-STATE-SHELLS-218
## The Mori-Zwanzig correction suppresses raw boundary-flux rates by factors of roughly four to ten

Date: 2026-09-26

Status:
exact finite-state M0/M1 calculation for N=6 and N=8.

Invert the instantaneous projected Fourier operator PLP into shell rates
q_d^(PLP), then repeat after memory renormalization.

At N=6 the ratios

    q_d^(MZ)/q_d^(PLP)

for d=1,...,6 are approximately:


    d=1: 0.094179

    d=2: 0.182996

    d=3: 0.234185

    d=4: 0.244537

    d=5: 0.189243

    d=6: 0.196333


At N=8:


    d=1: 0.117107

    d=2: 0.155418

    d=3: 0.194127

    d=4: 0.210159

    d=5: 0.161803

    d=6: 0.174248


Thus PLP overestimates persistent localized-to-localized motion severely.

The zero-frequency memory self-energy removes most of that raw flux before the
small first-moment clock correction is applied.

This gives a shell-resolved interpretation of the earlier recrossing result:

    PLP
      ≈ immediate boundary/reactive flux,

    M0
      ≈ return/recrossing self-energy,

    M1
      ≈ finite memory-time correction,

    Q_eff
      ≈ committed long-time localized-state motion.

The suppression is not uniform across shells, so memory changes not only the
overall clock but also the relative transition geometry.
