# A bounded corollary: a vanishing signal needs a paid readout error

Date:2026-10-15. Analytical consequence of the accepted SRN native bounds.
This is a diagnostic observation class, not an added FIN microscopic law
or an assertion that a real apparatus has this error.

Suppose the two reported channels are

    O_ri=a_i(W_ri+e_ri), r=1,3,

with a_i finite and nonzero, common to both channels at sample i. On each
uncorrupted sample suppose |e_ri|<=delta_N. Allow arbitrary dependence
on the path. If instead errors are specified outside a_i, their bounds must
be divided by |a_i|. Nonzero gains alone give no uniform bound on this
normalization, since |a_i| may approach zero.

Since absolute value is1-Lipschitz, a changed raw angular decision requires

    ||W3|-|W1||<=2delta_N.

For real x,y, ||x|-|y||=min(|x-y|,|x+y|). The native anti-concentration
therefore gives, at each stationary sample,

    P(H_observed!=H_true)<=C(delta_N+k^(-1/2)).

No independence of errors or repeated observations is needed. If at most
a deterministic fraction r_N of the samples are arbitrarily corrupted or
erased, then

    |average H_observed-average H_true|
       <=r_N+average 1{uncorrupted raw event changed}.

Markov's inequality at a fixed fraction of the ideal gap yields an added
classification-error allowance

    C(delta_N+k^(-1/2))/alpha_N,

provided r_N is smaller than that gap fraction. In particular a sufficient
condition for the same consistency conclusion is

    delta_N=o(alpha_N), r_N=o(alpha_N).

Transfer the combined bounded final event to actual iid parents once as
in the SRN proof. Do not transfer an unbounded error norm by unspecified TV.
This result assumes pointwise or jointly licensed error envelopes; estimated
hardware bounds require their own confidence, dependence and custody proof.

A fixed error allowance does not vanish in this sufficient bound when
alpha_N->0. This is **failure of this guarantee**, not proof that all fixed
noise models or all tests are impossible. Independent gains are not covered
by calling them additive errors without bounding their state-dependent
normalized effect. Physical sample clocks, feature assignment and common
block identity also remain supplied interface premises.

The corollary explains why another ideal classifier cannot substitute for
an actual model-to-system bridge. It is not a new empirical result, a
necessary robustness threshold, or permission to manufacture calibration.
