#!/usr/bin/env python3
"""P2717/S1667: concrete pseudoscalar/chiral source candidate matrix.

P2716 identified the correct representation type for a future sign source: an
inversion-odd pseudoscalar/chiral signed datum with a nonzero value exported by
strict artifacts.  P2717 performs the next honest theoretical-physics audit: it
checks concrete pseudoscalar/chiral source classes that a field theorist would
try first, and refuses to promote any class unless all source criteria are met.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.json"
MD = GEN / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2716_PSEUDOSCALAR_ACCEPTANCE": GEN / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.json",
    "P2715_SCALAR_SOURCE_NO_GO": GEN / "p2715_s1665_aut_equivariant_scalar_source_to_orientation_torsor_no_go.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2687_ANISOTROPIC_SOURCE": GEN / "p2687_s1637_anisotropic_source_class_audit.json",
    "P2709_RELEASE_BACKSCAN": GEN / "p2709_s1659_release_8_1_to_9_3_breakthrough_unlock_backscan.json",
}

ACCEPTANCE_CRITERIA = [
    "inversion_odd_representation",
    "strict_artifact_exports_source_law",
    "nonzero_signed_value_exported",
    "sign_not_orientation_convention",
    "coupling_to_p2708_p2714_orientation_torsor_exported",
]

CANDIDATES = [
    {
        "candidate_id": "levi_civita_volume_orientation_density",
        "physics_type": "orientation pseudoscalar density",
        "inversion_odd_representation": True,
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": False,
        "sign_not_orientation_convention": False,
        "coupling_to_p2708_p2714_orientation_torsor_exported": False,
        "blocker": "A Levi-Civita density is sign-odd, but its sign is an orientation convention unless a strict law chooses the orientation.",
    },
    {
        "candidate_id": "pontryagin_or_chiral_anomaly_density",
        "physics_type": "topological/chiral pseudoscalar density",
        "inversion_odd_representation": True,
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": False,
        "sign_not_orientation_convention": True,
        "coupling_to_p2708_p2714_orientation_torsor_exported": False,
        "blocker": "A Pontryagin/anomaly-type density would be the right parity class, but no strict artifact exports such a nonzero signed density or its coupling to the Z12 torsor.",
    },
    {
        "candidate_id": "eta_or_spectral_asymmetry_index",
        "physics_type": "spectral/chiral asymmetry invariant",
        "inversion_odd_representation": True,
        "strict_artifact_exports_source_law": False,
        "nonzero_signed_value_exported": False,
        "sign_not_orientation_convention": True,
        "coupling_to_p2708_p2714_orientation_torsor_exported": False,
        "blocker": "A spectral asymmetry could source chirality, but no current strict artifact exports an eta/index sign coupled to the boundary-cocycle orientation torsor.",
    },
    {
        "candidate_id": "oriented_z12_cycle_cup_product",
        "physics_type": "discrete cohomological pseudoscalar candidate",
        "inversion_odd_representation": True,
        "strict_artifact_exports_source_law": True,
        "nonzero_signed_value_exported": False,
        "sign_not_orientation_convention": False,
        "coupling_to_p2708_p2714_orientation_torsor_exported": True,
        "blocker": "The oriented Z12 cycle is the existing P2708/P2714 torsor itself; it supplies the two signs but not a non-premise rule selecting one.",
    },
]

NEGATIVE_EXPORT_FLAGS = [
    "concrete_pseudoscalar_source_accepted",
    "strict_artifact_exports_source_law",
    "nonzero_signed_value_exported",
    "sign_not_orientation_convention",
    "coupling_to_orientation_torsor_exported_as_closure",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def evaluate_candidates() -> list[dict[str, Any]]:
    rows = []
    for candidate in CANDIDATES:
        missing = [criterion for criterion in ACCEPTANCE_CRITERIA if not candidate[criterion]]
        accepted = not missing
        rows.append({
            **candidate,
            "missing_criteria": missing,
            "accepted_as_strict_source": accepted,
            "fixes_lambda_now": accepted,
        })
    return rows


def prior_statuses() -> dict[str, Any]:
    return {key: read_json(path).get("status") for key, path in INPUTS.items()}


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2717/S1667 concrete pseudoscalar/chiral source candidate matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate matrix",
    ]
    for row in payload["candidate_matrix"]:
        lines.append(f"- `{row['candidate_id']}`: accepted={row['accepted_as_strict_source']}, missing={row['missing_criteria']}. {row['blocker']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    rows = evaluate_candidates()
    accepted = [row for row in rows if row["accepted_as_strict_source"]]
    no_unlock = len(rows) == 4 and not accepted
    payload = {
        "status": "P2717_CONCRETE_PSEUDOSCALAR_CHIRAL_SOURCE_MATRIX_NO_ACCEPTED_SOURCE" if no_unlock else "P2717_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": prior_statuses(),
        "acceptance_criteria": ACCEPTANCE_CRITERIA,
        "candidate_matrix": rows,
        "accepted_candidate_count": len(accepted),
        "decision": {
            "concrete_pseudoscalar_source_accepted": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2717 audits concrete pseudoscalar/chiral source classes after P2716.  Levi-Civita orientation density, Pontryagin/anomaly density, eta/spectral asymmetry, and the oriented Z12 cup-product candidate all have the right sign intuition or parity, but none satisfies all strict criteria: exported source law, nonzero signed value, nonconventional sign, and coupling to the P2708/P2714 orientation torsor.  Therefore no candidate fixes lambda or discharges QW-2191.",
            "next_honest_step": "Do not keep enumerating generic pseudoscalar names.  A next admissible move must provide one explicit formula/artifact computing a nonzero signed pseudoscalar/chiral value and its coupling to the orientation torsor, or pivot to a different genuinely new typed object outside the closed lanes.  Otherwise preserve the P2697-P2717 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2717/S1667 concrete pseudoscalar/chiral source candidate matrix", "## P2717/S1667 concrete pseudoscalar/chiral source candidate matrix\n\n`P2717/S1667` audits concrete pseudoscalar/chiral source classes after P2716.  Levi-Civita orientation density, Pontryagin/anomaly density, eta/spectral asymmetry, and an oriented Z12 cup-product candidate are tested against the strict source criteria: inversion-odd representation, exported source law, nonzero signed value, nonconventional sign, and coupling to the P2708/P2714 orientation torsor.  No candidate satisfies all criteria, so no strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2717/S1667 pseudoscalar candidate matrix Ltotal guard", "## P2717/S1667 pseudoscalar candidate matrix Ltotal guard\n\n`P2717/S1667` is a source-candidate obstruction matrix, not a variational source construction.  No concrete pseudoscalar/chiral candidate exports a nonzero signed value coupled to the orientation torsor, so this does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current concrete pseudoscalar/chiral source matrix guardrail (P2717/S1667, 2026-06-14)", "## Current concrete pseudoscalar/chiral source matrix guardrail (P2717/S1667, 2026-06-14)\n\n- P2717 audits concrete pseudoscalar/chiral source classes: Levi-Civita orientation density, Pontryagin/anomaly density, eta/spectral asymmetry, and oriented Z12 cup-product.\n- None satisfies all strict criteria: exported source law, nonzero signed value, nonconventional sign, and coupling to the P2708/P2714 orientation torsor; therefore no strict `lambda` fixing or `QW-2191` discharge is exported.\n- Do not continue enumerating generic pseudoscalar names as closure evidence; a next admissible move must provide one explicit formula/artifact computing a nonzero signed pseudoscalar/chiral value and its torsor coupling, or pivot to a different genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    main()
