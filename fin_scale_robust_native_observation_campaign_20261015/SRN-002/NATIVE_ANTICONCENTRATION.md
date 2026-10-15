# SRN-002 — A uniform *native* conditional anti-concentration inequality

**Statement.** Assume the inherited positive MT cylinder supplies n_j>=cN for each j=0,...,11 and the entire comparison support (same c>0), with 1<=k<=N/2, N in 12N. Conditional on the complete population counts n, the fixed block is a uniform k-subset of N distinct labelled slots (multivariate hypergeometric). For h_±(j)=f3(j)±f1(j), H_±=k^-1/2 sum_{a in B}h_±(xi_a), for all real z and eta>=0,

 sup_{n in cylinder,z} P(|H_±-z|<=eta | n) <= C_c (eta+k^-1/2).

Consequently P(dist((W1,W3), {|w3-w1|=0} union {|w3+w1|=0})<=C eta | n) <=C_c(eta+k^-1/2). All lattice atoms, including ties, are included. There is no quantitative CLT premise.

**Exact two-type choices.** For + take labels (0,6), with h_+(0)=2sqrt(2),h_+(6)=-2sqrt(2), spacing 4sqrt(2). For − take (2,8), with h_-(2)=-3sqrt(2)/2,h_-(8)=+3sqrt(2)/2, spacing 3sqrt(2). Each type has at least cN original slots and the union has H=a+b>=2cN.

**Good-event mass.** Let t be the total number of chosen block slots from this two-type union. Conditional on n, t~Hypergeom(N,H,k), mean kH/N and variance <= k/4. Because H/N>=2c and k<=N/2, Chebyshev yields
 P(t<c k) <= C_c/k,   P(t>3H/4) <= C_c/k.
For the second bound E t<=H/2 and H>=2cN>=4ck, so (3H/4-Et)>=H/4>=c k. Thus with probability >=1-C_c/k, t>=ck and t/H<=3/4. For bounded k the claimed C_c k^-1/2 term absorbs trivial probabilities.

**Conditional one-dimensional lattice.** Conditional on t and every nonselected block count, the number X of first-type labels is hypergeometric (H,a,t):
 P(X=x|t,others,n)=[C(a,x)C(b,t-x)]/C(H,t).
The remaining labels contribute a deterministic shift to the sum; the selected pair contributes t h(second)+X(h(first)-h(second)). Consecutive values of the normalized total differ by Delta/sqrt(k), with |Delta|=4sqrt(2) or 3sqrt(2). Thus any interval |H_±-z|<=eta contains no more than 1+2eta sqrt(k)/|Delta| possible X-values (up to ceiling constants).

**Uniform maximum atom.** On the good event a/H,b/H>=c and t/H<=3/4 with t>=ck, the mode x_* of the above hypergeometric lies within O(1) of t a/H (use the ratio P(x+1)/P(x)=[(a-x)(t-x)]/[(x+1)(b-t+x+1)]). Each of x_*, t-x_*, a-x_*, b-t+x_* is >=c'_c t for t beyond a fixed threshold; for example a-x_*>=a(1-t/H)-O(1)>=cH/4-O(1). Apply the two-sided factorial/Stirling inequalities to the three binomial coefficients. The entropy/exponential factor at the mode is bounded by exp(C_c), since the mode differs O(1) from the entropy-maximizing table with proportions a/H and t/H; the prefactor is at most
 C_c sqrt(H^3/[t(H-t)ab]) <= C_c/sqrt(t) <=C_c/sqrt(k).
For finitely many small t enlarge C_c. Unimodality implies sup_x P(X=x|t,others,n)<=C_c/sqrt(k).

Multiplying max-atom by the number of lattice points gives C_c(eta+k^-1/2), uniformly in z. Add bad-event C_c/k <=C_c/sqrt(k). A union bound pays both angular boundaries. This is conditional *native* anti-concentration, not a finite-N independence approximation or an imported Gaussian density bound.
