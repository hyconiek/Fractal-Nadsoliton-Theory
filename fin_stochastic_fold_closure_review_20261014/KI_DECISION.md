# Architect decision after the SF intake

Date: 2026-10-14.
Decision: **KINETIC_INVARIANTS_AND_FOLD_RESPONSE**.
Status: next programme selected; KI execution has not started.

The prepared compact-stopped heat-bath bridge is sufficient for its next
theoretical layer after the intake's explicit joint-clock completion. No
extra static-root digits or another same-law simulation campaign is the
default. Numerical finite-N error constants and broader unprepared/global
limits remain interesting, but do not block this proposed diagnostic lane.

## The scientifically discriminating next question

> Which part of the critical trajectory law is fixed by the equilibrium FIN
> landscape, and which part depends on the supplied reversible kinetic rule?

FIN currently assumes its microscopic law. Earlier PHYS-004/E04 results show
that the same Gibbs equilibrium does not determine finite-time relaxation.
The new controlled critical law permits a much sharper question: can
thermodynamic coefficients and normalized response survive that kinetic
freedom, with mobility/noise derived rather than refitted?

Keep the primary FIN heat bath unchanged. Use only already specified
Metropolis/Barker rules as explicitly marked **countermodels**, with E04's
g0 per-copy clock convention. Extending their written local rate formulas to
an N-copy comparison family must be explicit and directly checked; N2 output
arrays are not a general-N proof. This introduces no new fundamental FIN law.

## A concrete route, still a proposed theorem

Use the primal occupation free energy F and the same exact fold right vector r.
At equilibrium, reversible edge conductances define a positive mobility K on
the normalized occupation tangent. The expected general relation is

    J=-K H_F,   D=2K,
    m_K=1/(r^T K_T^(-1) r),
    ell_K^T=m_K r^T K_T^(-1),
    dZ_K=-m_K(a lambda+b Z_K^2/2)ds+sqrt(2m_K)dB.

Here inverse means restriction to the tangent, not inverse of a singular
12x12 conserved-normalization matrix. For the primary heat bath K=Sigma_f,
m_K=g_f, recovering the admitted SF coefficients. The goal is to prove the
relation and limit for declared controls; **not** to promote this displayed
candidate as a completed cross-kinetic theorem.

There is an important new gate: the original thermal soft readout and the
kinetic left-null coordinate need not coincide. An integrated second-moment
estimate does not automatically make `N^(1/3)w_N` uniformly vanish along a
whole path. A fourth-moment generator estimate may pay that transfer for the
rounded prepared class. Metropolis also has equal-logit kinks at this fold;
smooth heat-bath expansions must not be copied silently.

These are useful failure tests, not reasons to invent an extra bath or fit a
clock to the desired canonical curve.

## Why this is preferable to another precision programme

- It tests how much physical content can survive unspecified kinetics.
- It distinguishes an equilibrium mechanism from a selected controller.
- It can prove universality or expose a surviving kinetic dependence.
- It prevents the normalized number 0.41685565... being sold as unique to A7.
- It reuses the completed source-to-fluctuation bridge without new experiments.

Begin KI-001, KI-002 and KI-003 from
[the English execution plan](../FIN_KINETIC_INVARIANTS_NEXT_CHAIN_20261014.md).
KI-004 is active only after its theorem dependencies pass; KI-005 returns the
handoff. Architect alone chooses a later source, geometry, global or empirical
branch. Missing source/locality/units/conservation gates remain unchanged.
