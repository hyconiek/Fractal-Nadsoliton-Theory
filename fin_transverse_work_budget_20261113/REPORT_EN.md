# Signed transverse reference-energy work budget

## Scope and identity

This repeats the conditional j=6 variational experiment at w=1.25,
epsilon=0.1, background psi phase pi/2, T=200 and periodic [-230,230).
The two fixed transverse seeds, evolving nonlinear background, reference
energy and all coefficients are unchanged from the horizon package.

For the fixed reference energy E_ref and z_tt=z_xx-V(t)z,
periodic integration by parts gives the exact identity

    dE_ref/dt = integral[(V_ref - V(t)) z z_t],
    V_ref - V(t) = -12 K eta - 6 eta^2 - 3 psi^2 - 5 psi^4,
    K=tanh(x).

Define P_phi=integral[(-12 K eta-6 eta^2) z z_t] and
P_psi=integral[(-3 psi^2-5 psi^4) z z_t]. These are signed powers in
the fixed reference norm. Even a nonpositive potential coefficient does
not determine the sign of power, because z z_t changes sign.

Both powers are integrated independently at the four RK4 stages, using
the same stage fields as the evolution. They are not reconstructed from
energy differences. The measured energy is calculated separately from
spatial derivatives and velocities. Every half-unit diagnostic stores
E_ref/E_initial and both cumulative works divided by E_initial for each
fixed seed. The residual is gain-1-W_phi-W_psi.

No work attribution is made to the time-varying optimizing Gram eigenvector:
that would mix changing initial conditions. Attribution uses the fixed
position and fixed velocity seeds throughout.

## Controls

Three runs use dx=0.1: static dt=0.02, driven dt=0.02 and driven dt=0.01.
Each energy trace is checked against the historical horizon run. Static
work must vanish. The analysis tests the balance at every stored time and
compares the complete cumulative-work traces between time steps. Spatial
refinement is inherited only from the earlier T=100 experiments; it has
not been repeated at T=200. No certified continuum error bound is claimed.
The detailed output is in comparison.json and RESULTS.md.

## Interpretation limits

This budget describes energy change in a chosen positive reference norm
for infinitesimal perturbations on a prescribed nonlinear background.
There is no back reaction, so these are not equal-and-opposite measured
losses of the background energy. They are not emitted radiation powers,
physical heating or separately conserved field energies. The split is
algebraically defined relative to the static kink and this field chart.

The result cannot establish nonlinear stability, an asymptotic decay law,
a strict FIN coefficient source, AW extension, units or empirical physics.
It explains the observed finite-window transfer within the supplied model.

## Replay

    OPENBLAS_NUM_THREADS=1 python fin_transverse_work_budget_20261113/run.py
    python fin_transverse_work_budget_20261113/analyze.py

Run from the repository root with Python and NumPy. Historical data are
read only; input and output hashes are retained.
