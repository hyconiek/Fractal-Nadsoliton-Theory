# Exact N2 observation process and a sharper design certificate

Date: 2026-10-14. Status: finite-model theorems and conditional interval
certificates added by the repository intake, not claimed as physical data.

## 1. Exact difference process

Assume the declared N2 parent with symmetric circulant A, zero field,
q=12, and the three specified label-equivariant kinetic laws. Write
a_d=A_(0,d), p_d=softmax_d[(g/2)a_d] and d=(sigma2-sigma1) mod12.
For every microscopic state in a fixed difference fibre, the sum of rates
to any other difference v is the same. For v!=d it is

| Supplied kinetics | Difference-generator rate B_(d,v) |
|---|---|
| Heat bath, one refresh per copy | 2 p_v |
| Metropolis, calibrated distinct-label proposal | (2/12) min(1,exp[(g/2)(a_v-a_d)]) |
| Barker, calibrated distinct-label proposal | (4/12)/(1+exp[-(g/2)(a_v-a_d)]) |

Choose diagonal entries to make rows sum to zero. This is the sum of the
two allowed single-copy moves to difference v. Circulant symmetry and
a_d=a_(-d) make their weights independent of absolute labels. All rates
are nonnegative and each law preserves p by detailed balance.

Consequently the partition is **strongly lumpable**, not an approximate
Galerkin closure. If R maps the 144 microscopic states to their difference,

    L R = R B,    exp(tL) R = R exp(tB).

This holds for every initial microscopic law. The three g=0 quotient
generators coincide, with off-diagonal rates 2/12. No time scale is fitted
to the g=3 outputs.

This theorem is separate from the nonlumpable M3 first-two-label observation.
It neither closes that process nor contradicts its explicit memory result.
d is an internal label difference, not a physical position or distance.

## 2. Heat-bath reset law

For heat bath,

    B = 2(1 p^T-I),
    exp(tB) = exp(-2t) I + (1-exp(-2t)) 1 p^T.

Thus all eleven nonconstant quotient modes decay at rate2, independent of
the particular p/A/g. This is a **kinetic-class** property, not FIN uniqueness.
The only stationary eigenvalue is0. For any observable f at equilibrium,

    Cov(f(d_t),f(d_0)) = exp(-2t) Var_p(f).

For an initial difference law mu0 and t2>t1,

    mu(t1) = exp(-2t1) mu0 + (1-exp(-2t1)) p,
    J_(d,v) = mu_d(t1) [exp(-2(t2-t1)) delta_(d,v)
                         +(1-exp(-2(t2-t1))) p_v].

In particular `TV(mu(t),p)=exp(-2t) TV(mu0,p)`. This pays mathematical
heat-bath burn-in; its conversion to a physical waiting time still needs
the calibrated production clock. The stationary separation certificate
cannot be applied directly to the source's early uniform-quench marginals.

The source's uniform quench uses mu0=1/12, t1=.4, t2=1. The closed formula
and all three 12-state generators reproduce all nine saved 13x13 ideal
joint laws, with the invalid row/column zero. Maximum numerical discrepancy
is less than 1.4e-16 in the intake. Matrix exponentials remain numerical;
the underlying algebra is exact.

For unequal constant copy-refresh rates the heat-bath factor2 is replaced
by their sum. Difference-only observations do not identify the two rates
separately. A physical clock and its invariance when gain is switched remain
independent calibration premises.

## 3. Observation instruments and scope limits

If a noninvasive emission channel depends only on d, with probabilities
k_o(d), the unnormalized column update is

    mu -> diag(k_o) exp(tB^T) mu.

Retaining all outputs, including INVALID, gives one positive, compositional
process for arbitrary record strings and arbitrary initial mu. The **noisy
output sequence itself need not be Markov**; the latent difference process
is Markov and its filter retains the required state information.

A general instrument on absolute labels, measurement disturbance,
noncirculant interaction, nonuniform label-dependent rates or fields does
not inherit this quotient. A symmetry-broken 0.01 off-diagonal perturbation
in the intake produces a nonzero within-fibre transition-rate defect. This
is a mathematical negative control, not a measured hardware defect.
Retain the 144-state parent or an enlarged physical state when needed.

Uniform g=0 histograms alone cannot certify readout: the constant emission
channel k_o(d)=1/12 for the twelve valid outputs reproduces the uniform
control while erasing every g=3 interaction fingerprint. Conditional/paired
calibration or a separately bounded observation channel is necessary.

## 4. The source's zero lower bound is not a separation no-go

The source used the triangle inequality and separate worst-case gain-TV
caps. Its FIN/flat lower bound is0 under hypothetical ±1% gain and .005
detector-TV allowance per model. That computation is valid but inconclusive:
it does not prove intersection of the underlying gain-constrained law sets.

The intake directly encloses **independently varying** gains

    g_FIN in [2.97,3.03],   g_FLAT in [2.97,3.03].

Each axis is covered by60 closed intervals of width .001. Exact-decimal
input and normalized probabilities are enclosed with mpmath.iv, with public
endpoints rounded outward by exact rational arithmetic. Every Cartesian
cell receives the lower bound

    TV >= (1/2) sum_d max(0, p_d_lower-q_d_upper,
                              q_d_lower-p_d_upper).

The complete **3,600-cell** certificate proves the conservative public result

    inf_(g_FIN,g_FLAT) TV(p_FIN,p_FLAT) >= 0.01119114.

A point witness at (2.97,2.97343) gives an upper bound about 0.011501225,
not an exact minimizing pair or exact minimax solution.

Embed the ideal law into thirteen outcomes with zero invalid mass. If each
allowed detector-output law is within TV .005 of its own ideal law, the
triangle inequality gives

    inf observed-law separation >= 0.00119114 > 0.

This allows independent gains and even independent output-error choices
within the stated balls. It does not rely on silently requiring the same
unknown detector parameter in both explanations. Nor does it include extra
unknown fields, interaction defects, drift or arbitrary detector errors.
Those must enter a larger explicit class and can remove this margin.

See [GAIN_BOX_CERTIFICATE.json](GAIN_BOX_CERTIFICATE.json). The standard-library
[checker](check_gain_certificate.py) recomputes the serialized coverage and
TV bounds in exact fractions and rejects an omitted coverage row.
Input probability boxes remain conditional on mpmath.iv correctness.

## 5. What the refinement permits

The stated hypothetical design class has a positive guaranteed separation.
The old loose-bound NOT_GUARANTEED status is superseded **for this class**.
The unconditional arbitrary-readout erasure no-go remains true. No actual
calibration has established that any physical apparatus lies in the class.

The proper next question is how to certify the instrument, symmetry,
preparation and clock premises operationally. Improving further digits of
this static margin is not the default next task.
