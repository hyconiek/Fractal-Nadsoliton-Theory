# Exact comparison seed: same parity mechanism, different fast spectrum

Date: 2026-10-15. Status: exact algebra under the inherited spectral boxes.
This defines **diagnostic comparison operators**, not replacement FIN laws.
Certificate: [KERNEL_FAMILY_SEED.json](KERNEL_FAMILY_SEED.json).
Checker: [kernel_family_certificate.py](kernel_family_certificate.py).

For Fourier projectors P3 and P4 on Z12, set

    A_epsilon=A7+epsilon(P3-P4),  |epsilon|<=1/100.

Both projectors have rank2 and diagonal1/6. Their entries are rational:
`P3(i,j)=cos(pi(i-j)/2)/6`, `P4(i,j)=cos(2pi(i-j)/3)/6`.
The checker evaluates their cosine tables in exact Fractions and verifies
idempotence, orthogonality, symmetry, centering and annihilation of parity.

Consequently the perturbation is nonzero but has zero diagonal/trace and
annihilates both uniform and parity directions. The family retains q12,
rank7, D12 covariance, positive active spectrum and the same unique top
eigenvalue Lambda6. Positivity and gap follow from exact rational comparisons
with the already accepted spectral boxes. No new transcendental point
calculation is promoted to an enclosure.

The changes are

    Lambda3 -> Lambda3+epsilon,
    Lambda4 -> Lambda4-epsilon,
    Lambda5,Lambda6 unchanged.

Thus g_c=12/Lambda6, trace, A00 and the **entire parity-line potential** are
identical, while fast curvatures differ. This is a particularly clean
comparison: no normalization or extra clock favors one operator.

What is not paid here: transfer of every MT/UF/LI theorem across this family,
uniform-in-epsilon limits, equality of full parity observation experiments,
or new fast readout laws. The KO plan pays only the structural conditions it
needs, first for epsilon0 and+1/100. The family does not derive nature's A7.

Critically, LI proves equality of risks for **one within-model future-sign
target**. Even if two operators have the same leading parity diffusion, this
does not prove that their complete microscopic parity histories are unable
to discriminate the operator. Shared weak limits give no such general no-go.
Only a fixed finite-resolution record can be declared asymptotically blind
after its full finite-alphabet distributions are proved to share a limit.

The useful next observable is already in FIN: fast Fourier occupation
fluctuations. MT's Gaussian fast precision contains `12I-g_c A`, so those
variables can retain spectral ratios erased from the canonical parity law.
Their static covariance is inherited mathematics, not a newly discovered
law. The new problem is a controlled partial readout and a bounded contrast.
