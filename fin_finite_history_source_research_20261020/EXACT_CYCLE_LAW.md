# FH-001: exact finite law, root convention and causal completions

Date: 2026-10-20. Conditional finite mathematics; known general method.

## 1. The canonical product is constant on each occupation class

The unchanged selected law satisfies, on legal edges,

    P_xy=s_x z_y/z_x.

For any legal L-step path with departure counts n,

    product_(t<L) P_(x_t,x_(t+1))
      =(z_(x_L)/z_(x_0)) product_x s_x^(n_x).          (1)

For a closed path the boundary ratio is one. Therefore the **trace-loop**
law proportional to product P, conditioned on n, is exactly U_n. It is
not the ordinary stationary path law conditioned on closing and n: that
law includes an additional pi_(x_0), so it biases the root. With a fixed
root a this factor is constant, and the two conditional laws agree.

Cyclic rotation is a bijection of the rooted word set. Averaging the exact
occupation identity over rotations proves

    U_n(X_t=x)=n_x/L                                  (2)

for every t, including periodic words. This is finite circular stationarity,
not independence or Markovity. Reversal is another such bijection; this
construction supplies no preferred arrow of time.

## 2. Exact next-state law without future guessing

Fix root a. Let C_a(x,r) be the number of completions from current x to a
with r_x departures from each state still to be consumed. For integer
r>=0 set

    C_a(x,0)=1{x=a},
    C_a(x,r)=0 if sum r>0 and r_x=0,
    C_a(x,r)=sum_y K_xy C_a(y,r-e_x) otherwise.        (3)

The recursion terminates since sum r drops by one. It counts legal paths
by their first step, so is exact, including self transitions. If C>0,

    Pr(next=y | x,r,a)=K_xy C_a(y,r-e_x)/C_a(x,r).     (4)

These probabilities sum to one. Their product telescopes to 1/C_a(a,n),
yielding the uniform rooted cycle. Choose a initially with probability
C_a(a,n)/sum_b C_b(b,n) to obtain U_n. By (2) this is n_a/L whenever
the class is nonempty. Zero completions are an impossible event, not an
invalid conditional distribution to evaluate or silently normalize.

The current state alone generally has memory. On the complete two-state
graph take root 0, n=(3,2) and L=5. At time 3 the two possible prefixes
(0,0,1,0) and (0,1,1,0) leave counts (1,1) and (2,0), respectively. At
the same current state 0 the former forces next=1, the latter next=0.
The variable r is an analysis record
of a finite conditioned ensemble, not a postulated physical controller,
new hidden substrate or additional observable supplied by FIN.

## 3. Root bias and scope controls

On a complete graph the selected P has every row pi. Given n, all closed
transition products are product pi_x^(n_x), but a stationary Markov closing
law weights rooted words by pi_(x_0). Its root marginal is proportional to
n_x pi_x, generally not n_x. Dropping that factor is not allowed.

For weighted K, equation (1) includes product K_(x_t,x_(t+1)); plain uniform
counting is then wrong. In particular the existing count-orbit reference
has off-diagonal sqrt(n_a(n_b+1)). Uniform labelled histories must be
projected with their actual multiplicities, not replaced by uniform count
words. No attempt marks are added to the one unmarked self continuation.
