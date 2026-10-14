# TH-003 path contraction

By TH-002, p_N obeys a path LDP with rate `I[p(.)]`. For the continuous observation `mu=X^T p`, the contraction principle gives

`I_mu[m(.)]=inf { I[p(.)] : X^T p(t)=m(t) for all t }`.

The infimum is over **entire hidden paths and the compatible initial-law cost**. An instantaneous minimization of the hidden state at each time is not generally equivalent and is not used.

Because `H(p,X eta)` varies inside a fibre of fixed mu (the N=120 witness already separates second derivatives), there is no universal exact Markov Hamiltonian on mu alone for arbitrary preparation. Three scopes must remain separate:

1. finite-N exact nonclosure;
2. local/stationary Gaussian reductions where hidden fluctuations can be slaved or averaged;
3. higher-order/nonstationary history dependence.

If exact path prediction is required without a reduction theorem, use full p or the exact 11-coordinate lift from `OBSERVABLE_SUFFICIENCY.md`.
