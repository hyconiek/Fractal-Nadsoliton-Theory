# DIRECT-SADDLE-VS-COMMUNICATION-CORRECTION-225
## Effective shell rates must be compared with communication paths, not automatically with same-distance direct saddles

Date: 2026-09-26

Status:
correction to the simplistic interpretation in reports 206/215.

Earlier comparisons placed:
- q3 against B3;
- q4 against B4;
- q5 against B5.

The first two are structurally reasonable because d3 and d4 are the edges that
create the two nontrivial levels of the barrier filtration.

The q5 comparison is not.

## 1. Why d5 is different

Although a direct d5 index-one saddle exists at

    B5=0.782654709774,

the two minima can communicate through a sequence of d3/d4 saddles whose
maximum barrier is only

    B4=0.662219137127.

Thus the communication height is B4, not B5.

The same issue applies to effective non-nearest shell rates such as d1,d2,d6:
their long-time influence can be induced by multi-step lower-barrier paths.

## 2. Consequence for q_d

The q_d reconstructed from low-frequency eigenvalues are EFFECTIVE Markov
rates.

They already contain:
- unresolved transition-layer elimination;
- recrossing renormalization;
- possible compression of multi-step communication.

Therefore q_d should not be assigned one-to-one to a direct same-distance
saddle unless that identification is independently proved.

## 3. Correct barrier statement

The robust current hierarchy is:

    intra-mod3 communication:
      B3;

    inter-mod3 communication:
      B4.

The direct B5 saddle is a higher alternative channel, not the bottleneck.

So the physically meaningful large-N candidate is a TWO-LEVEL communication
hierarchy, not six independent Eyring-Kramers shell laws.

## 4. Updated interpretation of finite-N q-rates

The observations

    q3 > q4 > q5

remain true for N=3..8.

But only the broad distinction

    mod3-internal vs mod3-external

should currently be tied to the B3/B4 barrier split.

Shell-by-shell exponent matching is deferred until a transition-path /
capacity theorem isolates genuinely direct channels.
