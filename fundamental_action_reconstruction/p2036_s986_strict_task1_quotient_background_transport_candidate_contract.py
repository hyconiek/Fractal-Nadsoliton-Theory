#!/usr/bin/env python3
"""P2036 S986: strict Task-1 quotient background-transport candidate contract.

P2035 proves that the P2034 local scalar B1 quotient class is not currently
transportable to background-global renormalization.  P2036 records the
minimal contract a future transport theorem would have to satisfy.

This packet deliberately exports a candidate contract, not a theorem.  It is a
new-premise work item until the required maps, finite-part scheme transport,
and global cocycle witnesses are actually computed.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.json"
MD = GEN / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.md"

SCHEMA_VERSION = "p2036_s986_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
NULL_DIRECTION = [1.0, -4.0, 1.0, -1.0]
QUOTIENT_BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]
Q_MATRIX_B1 = [
    [1.0, 0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0, -4.0],
    [0.0, 0.0, 1.0, 1.0],
]
IDENTITY_3 = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(value: Any) -> bool:
    return bool(value is True)


def det3(matrix: list[list[float]]) -> float:
    a = matrix
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def gate(gate_id: str, requirement: str, current_status: str, blocking_source: str) -> dict[str, str]:
    return {
        "gate_id": gate_id,
        "requirement": requirement,
        "current_status": current_status,
        "blocking_source": blocking_source,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p1879 = load("p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json")
    p1935 = load("p1935_s885_strict_po3_epsilon_bound_and_scheme_transfer_attempt_probe.json")
    p2034 = load("p2034_s984_strict_task1_quotient_only_renormalization_theorem.json")
    p2035 = load("p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json")

    p2034_checks = p2034.get("gatekeeper_checks") or {}
    p2035_checks = p2035.get("gatekeeper_checks") or {}
    p2034_theorem = p2034.get("quotient_renormalization_theorem") or {}
    target_q = p2034_theorem.get("target_quotient_coefficients_R2bar_Ric2bar_Riem2bar") or []
    p2035_required = (p2035.get("required_transport_contract") or {}).get("required_missing_exports") or []

    p2034_local_source_ready = (
        p2034.get("local_verdict") == "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE"
        and as_bool(p2034_checks.get("renormalization_licensed_only_in_quotient"))
    )
    p2035_obstruction_ready = (
        p2035.get("obstruction_verdict") == "PASS_CURRENT_EXPORT_NONTRANSPORTABILITY_WITH_TRACE"
        and as_bool(p2035_checks.get("obstruction_theorem_pass"))
        and not as_bool(p2035_checks.get("background_globalization_exported"))
    )
    p1879_has_anisotropy_contract = "bianchiI_transport_ansatz" in p1879
    p1935_scheme_transfer_is_partial = (
        (p1935.get("scheme_independence_transfer_lemma_v1") or {}).get("status") == "ATTEMPT_PARTIAL"
    )

    acceptance_gates = [
        gate(
            "C1_QUOTIENT_SOURCE",
            "Use only P2034 quotient data with coefficient order R2_bar/Ric2_bar/Riem2_bar.",
            "SATISFIED_BY_P2034_LOCAL_SOURCE",
            "none",
        ),
        gate(
            "C2_BASIS_PRESERVING_MAP",
            "Export rank-3 maps M_bg<-B1 on the same quotient basis for FRW and BianchiI.",
            "OPEN_SYMBOLIC_PLACEHOLDER",
            "P2035 basis_preserving_quotient_map_between_backgrounds",
        ),
        gate(
            "C3_FINITE_PART_SCHEME_TRANSPORT",
            "Prove the finite-part scheme transport map preserves the operator basis and counterterm subtraction.",
            "OPEN_SYMBOLIC_PLACEHOLDER",
            "P1935 finite-part transport compatibility blocker",
        ),
        gate(
            "C4_GB_NULL_TRANSPORT",
            "Classify the GB null direction as transported topological data or prove it remains quotiented out.",
            "OPEN_TOPOLOGICAL_CLASSIFICATION_MISSING",
            "P2035 GB_null_direction_transport_or_topological_classification",
        ),
        gate(
            "C5_ANISOTROPIC_WITNESS_DATA",
            "Compute Bianchi-I one-loop divergence data fixing the sigma2 transport correction.",
            "OPEN_LOOP_DATA_MISSING",
            "P1879 anisotropic one-loop integrals missing",
        ),
        gate(
            "C6_GLOBAL_COCYCLE",
            "Prove M_kj M_ji = M_ki on triple overlaps for renormalized quotient data.",
            "OPEN_COCYCLE_THEOREM_MISSING",
            "P1678/P2035 global atlas cocycle blocker",
        ),
    ]

    all_acceptance_gates_passed = all(row["current_status"].startswith("SATISFIED") for row in acceptance_gates)
    syntactically_complete_contract = (
        p2034_local_source_ready
        and p2035_obstruction_ready
        and p1879_has_anisotropy_contract
        and p1935_scheme_transfer_is_partial
        and det3(IDENTITY_3) == 1.0
        and len(p2035_required) >= 6
    )

    candidate_contract = {
        "contract_id": "Task1_B1_quotient_background_transport_candidate_v1",
        "contract_status": "CANDIDATE_EXPORTED_NOT_A_THEOREM",
        "new_premise_status": "NEW_PREMISE_CANDIDATE_NOT_ADMITTED_STRICT_PROOF",
        "source_object": "P2034 local scalar B1 quotient class [a]_B1",
        "source_target_quotient_coefficients_R2bar_Ric2bar_Riem2bar": target_q,
        "quotient_basis": QUOTIENT_BASIS,
        "source_null_direction_R2_Ric2_Riem2_GB": NULL_DIRECTION,
        "quotient_projection_matrix_Q_B1_from_R4_to_R3": Q_MATRIX_B1,
        "zeroth_order_transport_seed": {
            "M_FRW_from_B1_at_sigma2_0": IDENTITY_3,
            "M_BianchiI_from_B1_at_sigma2_0": IDENTITY_3,
            "identity_seed_rank": 3,
            "identity_seed_determinant": det3(IDENTITY_3),
            "warning": "Identity at sigma2=0 is a boundary seed, not a transport theorem.",
        },
        "unfixed_transport_ansatz": {
            "M_FRW_from_B1": "I3 + finite_part_scheme_correction_F_FRW",
            "M_BianchiI_from_B1": "I3 + sigma2*A_BI + finite_part_scheme_correction_F_BI + O(sigma2^2)",
            "M_global_chart_j_from_i": "I3 + C_ji on quotient basis, subject to cocycle C constraints",
            "unknowns": [
                "A_BI_ij anisotropic one-loop correction",
                "F_FRW_ij finite-part scheme correction",
                "F_BI_ij finite-part scheme correction",
                "C_ji atlas overlap correction",
            ],
        },
        "admission_rule": (
            "The candidate becomes a strict transport theorem only if every acceptance gate C1-C6 "
            "is supplied by computed witnesses and the transported divergence targets agree on the "
            "same quotient basis."
        ),
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2036",
        "stage_id": "S986",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_CANDIDATE_CONTRACT_WITH_TRACE",
        "result_kind": (
            "OPEN_NEW_PREMISE_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_EXPORTED__NO_GLOBALIZATION"
            if syntactically_complete_contract
            else "OPEN_INCOMPLETE_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE"
        ),
        "route": "strict_only_with_explicit_new_premise_candidate",
        "strict_lane_assumptions": [
            "strict_kernel_only",
            "no_legacy_transfer",
            "P2034_local_quotient_source_only",
            "candidate_contract_not_theorem",
            "no_background_globalization_claimed",
            "no_tensor_component_claims",
            "no_QW2191_selector_claim",
        ],
        "depends_on": {
            "p1879_present": p1879.get("_missing") is None,
            "p1935_present": p1935.get("_missing") is None,
            "p2034_present": p2034.get("_missing") is None,
            "p2035_present": p2035.get("_missing") is None,
            "p2034_local_source_ready": p2034_local_source_ready,
            "p2035_obstruction_ready": p2035_obstruction_ready,
        },
        "input_hashes": {
            "p1879_json_sha256": file_sha256(GEN / "p1879_s829_strict_frw_to_bianchiI_transport_contract_probe.json"),
            "p1935_json_sha256": file_sha256(GEN / "p1935_s885_strict_po3_epsilon_bound_and_scheme_transfer_attempt_probe.json"),
            "p2034_json_sha256": file_sha256(GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json"),
            "p2035_json_sha256": file_sha256(GEN / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json"),
        },
        "grep_audit_decision": {
            "p2036_already_done_before_this_packet": False,
            "closest_existing_artifacts": [
                "P1879: transport ansatz for amplitudes, not quotient counterterms",
                "P1935: partial scheme-transfer attempt, not a contract for P2034 quotient data",
                "P2035: obstruction theorem listing required missing transport exports",
            ],
            "missing_before_this_packet": (
                "No explicit P2036/S986 candidate contract specified the quotient-basis maps and "
                "acceptance gates needed to upgrade P2034 beyond local B1."
            ),
        },
        "professor_decision": {
            "decision": "EXPORT_NONSTRICT_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_KEEP_TASK1_LOCAL",
            "rationale": [
                "P2035 blocks globalization from current exports but also names the missing map family.",
                "The productive next move is to make the future transport theorem checkable as a contract.",
                "The contract is explicitly a new-premise candidate until loop data, scheme transport, and cocycle witnesses are supplied.",
            ],
            "rejected_branch_for_now": "declare background-global renormalization from the identity seed alone",
        },
        "candidate_contract": candidate_contract,
        "acceptance_gates": acceptance_gates,
        "gatekeeper_checks": {
            "p2034_local_source_ready": p2034_local_source_ready,
            "p2035_obstruction_ready": p2035_obstruction_ready,
            "p1879_anisotropy_contract_available": p1879_has_anisotropy_contract,
            "p1935_scheme_transfer_still_partial": p1935_scheme_transfer_is_partial,
            "candidate_contract_syntactically_complete": syntactically_complete_contract,
            "identity_seed_rank3": det3(IDENTITY_3) == 1.0,
            "all_acceptance_gates_passed": all_acceptance_gates_passed,
            "strict_transport_theorem_exported": False,
            "background_globalization_exported": False,
            "finite_part_scheme_transport_proven": False,
            "gb_null_transport_classified": False,
            "global_cocycle_proven": False,
            "no_global_renormalization_claimed": True,
            "no_tensor_projection_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "theorem_scope": {
            "licensed_statement": (
                "P2036 licenses only a syntactic candidate contract for future background transport "
                "of the P2034 quotient data.  It does not license transport itself."
            ),
            "not_licensed": [
                "background-global Task-1 renormalization",
                "FRW or BianchiI transported counterterm values",
                "finite-part scheme transport theorem",
                "GB null direction topological transport theorem",
                "global atlas cocycle theorem",
                "tensor-component B1 renormalization",
                "QW-2191 selector closure",
                "ToE closure",
            ],
        },
        "false_pass_guard": (
            "The identity seed and symbolic matrices are contract scaffolding only.  They must not be "
            "used as evidence of background independence or global renormalization."
        ),
        "next_honest_step": (
            "Build P2037 by computing or importing one concrete same-scheme finite-part map on the "
            "R2_bar/Ric2_bar/Riem2_bar quotient basis, then test C2-C3 before touching global cocycles."
        ),
        "lay_explanation": (
            "P2036 robi formularz kontrolny dla przyszlego przenoszenia wyniku B1.  Wpisalismy, "
            "jakie mapy i testy beda potrzebne, ale pola sa nadal puste obliczeniowo, wiec wynik "
            "pozostaje lokalny."
        ),
        "environment": {
            "python": platform.python_version(),
        },
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    gate_lines = "\n".join(
        f"- `{row['gate_id']}`: `{row['current_status']}`" for row in acceptance_gates
    )
    not_licensed_lines = "\n".join(
        f"- {item}" for item in payload["theorem_scope"]["not_licensed"]
    )

    md = f"""# P2036 S986 Strict Task-1 Quotient Background Transport Candidate Contract

Status: `OPEN_CANDIDATE_CONTRACT_WITH_TRACE`

Result kind: `{payload['result_kind']}`

## Professor Decision

`EXPORT_NONSTRICT_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_KEEP_TASK1_LOCAL`

P2036 exports the contract shape needed after P2035.  It is explicitly a
new-premise candidate, not a strict transport theorem.

## Candidate Contract

Contract id: `Task1_B1_quotient_background_transport_candidate_v1`

Source object: `P2034 local scalar B1 quotient class [a]_B1`

Quotient basis:

`{QUOTIENT_BASIS}`

Null direction:

`(1, -4, 1, -1)`

Zeroth-order seed:

`M_FRW_from_B1 = I3`, `M_BianchiI_from_B1 = I3` at `sigma2=0`

This identity seed is only a boundary seed.  It is not background
independence.

## Acceptance Gates

{gate_lines}

## Not Licensed

{not_licensed_lines}

## Next Honest Step

Compute or import one concrete same-scheme finite-part map on the
`R2_bar/Ric2_bar/Riem2_bar` quotient basis, then test C2-C3.
"""
    MD.write_text(md, encoding="utf-8")
    print(OUT)
    print(MD)


if __name__ == "__main__":
    main()
