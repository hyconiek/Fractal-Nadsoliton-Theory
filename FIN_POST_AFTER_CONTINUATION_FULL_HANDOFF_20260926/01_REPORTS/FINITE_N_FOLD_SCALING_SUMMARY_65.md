# FINITE-N-FOLD-SCALING-SUMMARY-65
## A single certified saddle-node now controls barrier, noise and time scaling

Date: 2026-09-26

The current local fold package now has a closed leading scaling structure.

For

    delta=g-g_f,

the accepted/derived quantities are:

    branch displacement
      ~ delta^(1/2);

    slow rate
      ~ delta^(1/2);

    local barrier
      ~ delta^(3/2);

    N-copy noise amplitude
      ~ N^(-1/2).

Balancing barrier and noise gives

    delta~N^(-2/3).

Inside that window:

    branch displacement ~ N^(-1/3),
    soft standard deviation ~ N^(-1/3),
    relaxation time ~ N^(1/3),
    N Delta V = O(1).

The leading stochastic normal form is

    dY
      =
      -g_f[a Delta+(b/2)Y^2]d tau
      +sqrt(2g_f)dW_tau.

All coefficients except the overall physical event-rate conversion are inherited
from already declared/certified FIN objects.

The remaining research problem is no longer the local scaling exponent.
It is the GLOBAL question:

    Is the R7P-031 local saddle the relevant communication saddle of the full
    finite-N landscape, or can another path escape with a lower V_g barrier?

That should be addressed by a global communication-height search/certificate,
not by further local Taylor expansion.
