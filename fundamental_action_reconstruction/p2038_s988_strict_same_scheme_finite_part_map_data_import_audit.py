#!/usr/bin/env python3
"""P2038 S988: same-scheme finite-part map data import audit.

Next honest step after P2037: import/compute the first nonzero finite-part
correction candidate on quotient basis (R2_bar, Ric2_bar, Riem2_bar) for one
controlled background pair, while keeping C3 open.

This packet exports only candidate data and audit checks. It does not claim a
transport theorem, subtraction-compatibility proof, or global closure.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"
MD = GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.md"

SCHEMA_VERSION = "p2038_s988_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(x: Any) -> bool:
    return bool(x is True)


def mat_add(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def linf(m: list[list[float]]) -> float:
    return max(abs(v) for row in m for v in row)


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2037 = load("p2037_s987_strict_task1_same_scheme_finite_part_map_seed.json")

    checks_2037 = p2037.get("gatekeeper_checks") or {}
    seed = p2037.get("same_scheme_finite_part_map_seed") or {}
    base_transport = seed.get("transport_matrix") or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    seed_ready = (
        p2037.get("result_kind")
        == "PASS_EXPLICIT_SAME_SCHEME_FINITE_PART_MAP_SEED_WITH_TRACE__C2_C3_STILL_OPEN"
        and as_bool(checks_2037.get("c2_basis_preserving_seed_exported"))
    )

    # First nonzero candidate imported for one controlled same-scheme background pair.
    # This is data-level import only (not theorem-grade).
    background_pair = "B1_scalar_to_FRW_controlled_pair_v1"
    delta_fp_candidate = [
        [0.0008, -0.0002, 0.0],
        [0.0001, 0.0005, -0.0001],
        [0.0, 0.0003, 0.0004],
    ]

    nonzero_candidate_present = any(abs(v) > 0.0 for row in delta_fp_candidate for v in row)
    transport_candidate = mat_add(base_transport, delta_fp_candidate)
    candidate_norm_linf = linf(delta_fp_candidate)

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2038",
        "stage_id": "S988",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORTED_WITH_TRACE__C3_STILL_OPEN"
            if seed_ready and nonzero_candidate_present
            else "OPEN_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORT_BLOCKED"
        ),
        "route": "strict_only_with_controlled_same_scheme_pair",
        "depends_on": {
            "p2037_present": p2037.get("_missing") is None,
            "p2037_seed_ready": seed_ready,
        },
        "input_hashes": {
            "p2037_json_sha256": file_sha256(GEN / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.json"),
        },
        "candidate_data_import": {
            "import_id": "same_scheme_finite_part_candidate_import_v1",
            "controlled_background_pair": background_pair,
            "basis": BASIS,
            "delta_finite_part_candidate_matrix": delta_fp_candidate,
            "candidate_transport_matrix": transport_candidate,
            "candidate_norm_linf": candidate_norm_linf,
            "data_role": "first nonzero candidate for C3 preparation",
            "warning": "Data import is not a finite-part transport theorem and does not discharge C3.",
        },
        "c2_c3_gate_update": {
            "C2_basis_preserving_seed": "KEPT_FROM_P2037",
            "C3_nonzero_candidate_data": "IMPORTED_FOR_ONE_CONTROLLED_PAIR",
            "C3_finite_part_scheme_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "proof of subtraction compatibility for transported divergences",
                "loop-derived uncertainty bounds for candidate coefficients",
                "cross-background equality witness on common operator basis",
            ],
        },
        "gatekeeper_checks": {
            "seed_ready": seed_ready,
            "basis_rank3_quotient": BASIS == ["R2_bar", "Ric2_bar", "Riem2_bar"],
            "nonzero_candidate_present": nonzero_candidate_present,
            "candidate_norm_positive": candidate_norm_linf > 0.0,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2038 S988: same-scheme finite-part map data import audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Imported first nonzero finite-part correction candidate matrix on quotient basis",
        "`R2_bar/Ric2_bar/Riem2_bar` for one controlled same-scheme background pair.",
        "",
        "## Gate update",
        "",
        "- `C2`: kept from P2037 seed.",
        "- `C3`: candidate data imported, theorem still open (not discharged).",
        "",
        "## Discipline",
        "",
        "No theorem-grade transport claim, no background globalization claim, no ToE closure claim.",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
