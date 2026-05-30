# D12 chi_11 max shell-imbalance selector certificate

Status: `branch-generator-uniquely-maximizes-shell-labelled-d1-d5-imbalance-but-imports-axis`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- D12 orbit count: `38`
- Unit-5 two-cycle count: `13`

## Selector summary

- Candidate score: `abs(h_d5-h_d1) on each D12/unit5 two-cycle`
- Requires shell label: `True`
- Requires reduced D12 quotient: `True`
- Max |h_d5-h_d1|: `4`
- Maximizer count: `1`
- Unique maximizer is branch cycle: `True`
- Amplitude distribution: `{0: 2, 1: 6, 2: 4, 4: 1}`
- Full-Aut allowed strict source: `False`

## Branch maximizer

- Cycle index: `0`
- Orbit pair: `[0, 37]`
- Low/high histograms: `[4, 3, 2, 1, 0, 0]` / `[0, 3, 2, 1, 4, 0]`
- Low/high d5-d1: `-4` / `4`
- Low/high gap necklaces: `[1, 1, 1, 1, 8]` / `[2, 2, 3, 2, 3]`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are enumerated, quotiented by D12, and reduced to the unit-5 two-cycle basis.
- `covariance_fact`: The shell-labelled imbalance I=h_d5-h_d1 satisfies I(unit5*O)=-I(O) on every D12 two-cycle.
- `maximality_certificate`: Across the 13 chi_11 two-cycles, max |I| is 4 and is attained by 1 cycle(s); the amplitude distribution is {0: 2, 1: 6, 2: 4, 4: 1}.
- `branch_selection`: The unique maximizer is the branch cycle [0,37], i.e. A1/A11 versus A5/A7, so max |h_d5-h_d1| conditionally selects the branch generator.
- `import_boundary`: The rule is not full-Aut strict provenance because it uses both the D12 quotient and the labelled d1/d5 shell axis.
- `full_aut_intersection`: If full-Aut invariance is restored, d1 and d5 are exchanged and this oriented imbalance is not an invariant source.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result is conditional on D12 quotient plus shell-labelled d1/d5 data; it is not a full-Aut strict-source theorem.
