# P2018 S968 Strict Cutkosky P2017 Provenance Non-Promotion Audit Theorem

Status: `P2018_EXECUTED_P2017_P1953_CONTRACT_NONPROMOTION_AUDIT_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

Respond to the P2017 provenance concern without adding another physical ansatz.
P2018 audits the P2017 strict-kernel quadrature candidate against the earlier
P1953 dressed-amplitude interface contract and mechanically prevents P2017's
numerical diagnostic pass from being promoted into a backend-derived Cutkosky
theorem.

## Audited objects

P2018 audits:

```text
p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json
```

against the P1953 contract source:

```text
p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json
```

P1953 already specified the minimum theorem-grade dressed-amplitude interface:
common-basis dressed amplitude, squared amplitude, DiscM/CutSum, BRST projector,
dressed poles/residues, ghost trace, phase-space measure, and uncertainty or
exactness certificate.

## Non-promotion rule

The rule enforced by P2018 is:

```text
promote_to_backend_cutkosky_closure
  => (provenance_complete and p1953_contract_complete).
```

Equivalently, if any required provenance gate is false or any required P1953
contract field is missing/partial, theorem-level promotion is blocked even when
quadrature numerics pass.

## Machine-checkable audit

P2018 exports:

1. the P2017 result/status summary,
2. a field-by-field P1953 contract audit of P2017,
3. missing and partial P1953 contract fields,
4. the list of missing P2017 provenance gates,
5. the list of failed P2017 numerical gates, if any,
6. a numerical audit of the P2017 quadrature covariance candidate,
7. a symbolic non-promotion rule using SymPy implication,
8. gatekeeper checks proving that the promotion block is active.

The intended audit result is:

```text
PASS_PROVENANCE_NONPROMOTION_AUDIT
```

with global obstruction status still:

```text
OPEN_OBSTRUCTION_WITH_TRACE.
```

## No-false-pass boundary

P2018 does not prove unitarity.  It proves a narrower governance theorem:
P2017 must remain a quadrature-candidate diagnostic until it satisfies the P1953
same-scheme dressed-amplitude contract and exports strict-side amplitude
provenance.  It does not discharge `QW-2191`, prove all-state Cutkosky closure,
or close ToE.

## Next honest step

Build one strict-side dressed-amplitude interface component from the P1953
contract, starting with a same-scheme BRST-projected `M_dressed_common_basis`
component, then rerun P2018 to see which contract and provenance gates can be
truthfully upgraded.
