# MAXENT-BINARY-FIBER-REFRESH-124
## Maximum entropy fixes the binary update shape but leaves its event rate

Date: 2026-09-26

Status:
- exact finite-state theorem;
- direct binary specialization of reports 54-55.

Suppose a two-child fiber has stationary target

    q_f=(1/2,1/2).

Among all one-event kernels stationary at q_f, maximum conditional entropy
selects uniquely

    K_f
      =
      [[1/2,1/2],
       [1/2,1/2]].

So a fiber event completely forgets the old child and redraws ± uniformly.

If these events occur at Poisson rate

    rho_f,

the continuous-time generator is

    Q_f
      =
      rho_f(K_f-I)

      =
      -(rho_f/2)L_2.

Therefore the ST231 Laplacian coefficient is

    boxed:
    mu=rho_f/2.

And the odd fiber relaxation rate is

    2mu=rho_f.

Thus maximum entropy DOES remove the arbitrary acceptance/update-shape freedom,
but it does NOT select the Poisson event rate rho_f.

The refinement ambiguity has been reduced to a clock/activity ambiguity.

Three distinct source hypotheses remain possible:

1. independent fiber clock:
       rho_f is a new parameter;

2. inherited microscopic FIN refresh clock:
       rho_f=rho_micro;

3. inherited metastable slow clock:
       rho_f=delta_meta.

Only option 3 reproduces the conditional spectral-matching candidate

    mu=delta_meta/2.

None of these clock-identification rules follows from maximum entropy alone.
