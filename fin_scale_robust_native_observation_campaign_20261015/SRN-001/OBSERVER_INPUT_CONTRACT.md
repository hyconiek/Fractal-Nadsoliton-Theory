# SRN-001 — What the classifier actually sees

At i=0,...,M-1 the observer receives only a pair of block aggregates from the same predetermined tagged block: (W1,N(t_N+i),W3,N(t_N+i)), or their common-gain multiplied pair (a_i W1,a_i W3) with a_i!=0. It computes H_i=1{|W3|>|W1|} and A_N=M^-1 sum H_i. It does not receive total n, F1,F3,R1,R3, the q12 labels outside the block, the hidden parent model ID, the mollified test, Poisson corrector or the proof skeleton endpoint at time M.

Allowed nuisance: arbitrary finite nonzero common a_i for each i (even negative and varying). Not licensed: separate channel gains, offsets, erasure, channel noise, uncertain feature labels, missing observations or real-world instrument calibration. Neither a real detector nor a physical spatial half/block follows from mathematical tagging.
