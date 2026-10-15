# SRN-002 — Conditional factorial detail and conditioning scope

Given n=(n_0,...,n_11) and k, fixed-block assignments have probability prod_j C(n_j,b_j)/C(N,k). Select two types (a,b), with population sizes A=n_a, B=n_b and H=A+B. Condition on all outside-two-type block counts b_j, j not in {a,b}; their sum determines t=k-sum_{j other}b_j. In the remaining union, X=b_a has exactly the law C(A,X)C(B,t-X)/C(H,t). The conditional law does **not** pretend labels to be independent, and keeps all possible X including endpoints.

Ratio p(x+1)/p(x)=(A-x)(t-x)/[(x+1)(B-t+x+1)] is nonincreasing in x on its support; this gives unimodality and a mode within one of (t+1)(A+1)/(H+2) up to conventional integer rounding. In the good regime A/H,B/H>=c, t/H<=3/4, t>=ck, all four mode table cells grow at least c' t. With Robbins/Stirling's e^{1/(12n)} factorial remainder and binomial entropy expansion, the mass of the mode is bounded by C_c sqrt(H^3/[t(H-t)AB]). Since AB>=c^2H^2 and H-t>=H/4, this is <= C_c t^-1/2. The estimate is independent of N and of a shift caused by other labels.

This is a single fixed n conditional theorem; to obtain a stationary result average it over pi_D. It presupposes the inherited D_N positive-cylinder property. It does not transfer an L^p observable across a small-TV comparison or discard any boundary atom.
