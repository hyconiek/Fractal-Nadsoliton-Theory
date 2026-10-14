# Probability scope — RG-F01

**Exact finite-model result.** Let `R[x,d]=1{d=(b-a) mod12}`, for `x=(a,b)` on 144 labelled states. For every specified zero-field circulant, equivariant generator `L`, the rows satisfy `L R = R B`; `B` has nonnegative off-diagonal rates and zero row sums. Therefore a single 12-state CTMC works for **every microscopic initial phase distribution**, after pushing it to the difference marginal. The heat-bath `B=2(1 p^T-I)` and `exp(tB)=exp(-2t)I+(1-exp(-2t))1 p^T`; the other two B matrices are distinct. The process remains exact only for this difference-fibre model.

**Exact instrument rule.** A fixed column-stochastic 13×12 matrix `K[o,d]` defines a nondisturbing classical emission (including `INVALID`) with branch update `v -> v diag(K[o,:])` after propagation. Branches stay unnormalized. Summing over all outcomes gives the propagated law; summing intermediate observations recovers the shorter record law. This constructs one positive, compositional likelihood for arbitrary finite sequences and mixed initial laws. No basis/parameter is refitted per record.

**Physical assumption.** A detector that disturbs the state needs an independently calibrated family `J_o[d,d']` (or a full 144-state successor kernel). `K` is not sufficient to represent disturbance or an invalid physical state that evolves differently. An observed noisy output string need not itself be a Markov chain. Phase-dependent emissions may destroy the 12-state instrument representation even though the parent kinetics remain lumpable.

**Arithmetic.** Serialized matrices and laws are double-precision numerical witnesses; positivity of the abstract CTMC is an exact mathematical consequence of the nonnegative rate definition. No full-operation interval certificate is claimed.
