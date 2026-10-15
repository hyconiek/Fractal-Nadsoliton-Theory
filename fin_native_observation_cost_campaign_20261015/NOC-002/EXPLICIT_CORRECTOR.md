# NOC-002 — Frozen smooth reference and discrete Poisson corrector

Freeze `gamma in {gamma_A7,gamma_control}` with both strictly positive, `0<alpha<=1/2`, and `c=1-alpha` in the normalized private chart. The purely analytical OU operator is

`G=-gamma f partial_f + partial_ff -r partial_r +(1-alpha)partial_rr`.

Its time-t Gaussian transition has deterministic means `(exp(-gamma t) f, exp(-t) r)` and variances `(1-exp(-2gamma t))/gamma` and `(1-alpha)(1-exp(-2t))` independently. Therefore for `g_alpha(f,r)=cos(sqrt(alpha)f+r)`:

`P_t g_alpha(f,r) = exp(-[alpha(1-exp(-2gamma t))/gamma+(1-alpha)(1-exp(-2t))]/2) cos(sqrt(alpha)exp(-gamma t)f+exp(-t)r)`.

Stationary mean: `mu_gamma,alpha=exp(-[alpha/gamma+1-alpha]/2)`.

Let `psi_gamma,alpha(f,r)=sum_{ell=0}^infinity [P_ell g_alpha(f,r)-mu_gamma,alpha]`.

The series is absolutely convergent at every finite `(f,r)` and locally uniformly, and satisfies the exact identity `(I-P_1) psi = g_alpha-mu` by telescoping. This is a **discrete-time Poisson corrector**, for snapshot spacing exactly one model-time unit. It is an analysis function; the native generator remains the original jump generator.
