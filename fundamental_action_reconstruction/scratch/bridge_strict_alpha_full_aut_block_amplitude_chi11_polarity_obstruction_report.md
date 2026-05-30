# Full-Aut block amplitude and chi_11 polarity obstruction

Status: `full-aut-amplitude-uniquely-locates-branch-block-but-not-chi11-polarity`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- D12 orbit count: `38`
- Full affine Aut orbit count: `25`

## Block amplitude summary

- Candidate score: `max_component_abs(h_d5-h_d1) per full-Aut block`
- Full-Aut block score: `True`
- Exports chi_11 polarity: `False`
- Max amplitude: `4`
- Maximizer count: `1`
- Unique maximizer is branch full-Aut block: `True`
- Amplitude distribution: `{0: 14, 1: 6, 2: 4, 4: 1}`
- Branch block signed values: `[-4, 4]`

## Branch block certificate

- Full-Aut orbit index: `0`
- D12 components: `[0, 37]`
- Signed d5-d1 values: `[-4, 4]`
- Unoriented amplitude: `4`
- Has polarity pair: `True`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are enumerated and partitioned into D12 and full affine Aut orbits.
- `amplitude_invariance`: The quantity |h_d5-h_d1| is unchanged when unit 5 swaps d1 and d5, so it is a full-Aut block amplitude rather than a polarity.
- `unique_block_location`: Across the 25 full-Aut blocks, max amplitude is 4, attained by 1 block; amplitude distribution is {0: 14, 1: 6, 2: 4, 4: 1}.
- `polarity_obstruction`: Inside the unique branch block the two D12 components have signed values -4 and +4, so full-Aut gluing keeps only the block and erases the chi_11 sign.
- `selector_boundary`: The block can be located without choosing polarity, but A5 over A1 still requires a unit-axis or shell-orientation premise.
- `strict_limit`: This is not a strict-core chi_11 source theorem and does not discharge QW-2191.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result locates a full-Aut branch block but does not select chi_11 polarity or close the strict-core selector obstruction.
