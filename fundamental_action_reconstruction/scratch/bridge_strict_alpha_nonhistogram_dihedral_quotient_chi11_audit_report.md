# Non-histogram dihedral-quotient chi_11 audit

Status: `dihedral-nonhistogram-features-can-carry-chi11-only-with-reduced-symmetry-premise`

## Finite model

- Ring: `Z_12`
- Active support count: `5`
- Enumerated supports: `792`
- D12 orbit count: `38`
- Full affine Aut orbit count: `25`
- Full-Aut two-D12-component blocks: `13`
- Full-Aut fixed D12 blocks: `12`

## Branch certificate

- Branch full-Aut orbit index: `0`
- D12 component indices: `[0, 37]`
- Component gap necklaces: `[[1, 1, 1, 1, 8], [2, 2, 3, 2, 3]]`
- Unit 5 maps A1 component to A5 component: `True`
- Full-Aut invariant singleton A5-not-A1 classifier count: `0`

## Classifier counts

- D12 Boolean classifiers: `2^38 = 274877906944`
- D12 branch-chi_11-compatible classifiers: `2^36 = 68719476736`
- Full-Aut Boolean classifiers: `2^25 = 33554432`
- Full-Aut singleton A5-not-A1 classifiers: `0`
- Global integer chi_11 anti-invariant dimension on D12 quotient: `13`
- Global nonzero ±1 chi_11 character count: `0`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are enumerated exactly.
- `quotient_split`: The D12 quotient has 38 orbits; the full affine Aut quotient has 25 orbits.
- `branch_block`: A1/A11 and A5/A7 are two distinct D12 orbits but one full-Aut affine orbit; unit 5 swaps them.
- `full_aut_no_go`: Any full-Aut invariant support-only source is constant on the full affine orbit, so it has F(A1)=F(A5) and cannot export nonzero chi_11.
- `dihedral_positive_boundary`: A D12-invariant sign can be assigned as - on the A1/A11 component and + on the A5/A7 component, matching chi_11 on the branch block.
- `import_diagnosis`: That D12 sign is not strict full-Aut data: it requires a reduced cyclic/dihedral order, equivalently a unit-axis premise separating {±1} from {±5}.
- `global_character_space`: Under the unit-5 action on D12 orbits there are 13 two-component blocks and 12 fixed blocks; integer chi_11 anti-invariants have dimension 13 and must vanish on fixed blocks.
- `scope_limit`: This audits support-combinatorics features, not all possible strict nadsoliton provenance mechanisms.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result is a finite support-combinatorics boundary audit, not a complete strict-source no-go theorem.
