#!/usr/bin/env python3
"""P2037 S987: strict Task-1 first same-scheme finite-part map seed.

Implements the next honest move after P2036: export the first explicit
same-scheme finite-part map on the quotient basis (R2_bar, Ric2_bar, Riem2_bar).

This packet is intentionally minimal and local.  It exports a seed map for the
same renormalization scheme and records that C2/C3 are still not discharged:
no anisotropic loop witness and no theorem-grade finite-part transport proof.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.json"
MD = GEN / "p2037_s987_strict_task1_same_scheme_finite_part_map_seed.md"

SCHEMA_VERSION = "p2037_s987_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]
I3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]


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


def mat_sub(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[a[i][j] - b[i][j] for j in range(3)] for i in range(3)]


def linf(m: list[list[float]]) -> float:
    return max(abs(x) for row in m for x in row)


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2034 = load("p2034_s984_strict_task1_quotient_only_renormalization_theorem.json")
    p2035 = load("p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json")
    p2036 = load("p2036_s986_strict_task1_quotient_background_transport_candidate_contract.json")

    p2034_checks = p2034.get("gatekeeper_checks") or {}
    p2035_checks = p2035.get("gatekeeper_checks") or {}
    p2036_checks = p2036.get("gatekeeper_checks") or {}

    local_source_ready = (
        p2034.get("local_verdict") == "PASS_QUOTIENT_ONLY_RENORMALIZATION_WITH_TRACE"
        and as_bool(p2034_checks.get("renormalization_licensed_only_in_quotient"))
    )
    obstruction_context_ready = (
        p2035.get("obstruction_verdict") == "PASS_CURRENT_EXPORT_NONTRANSPORTABILITY_WITH_TRACE"
        and as_bool(p2035_checks.get("obstruction_theorem_pass"))
    )
    contract_context_ready = (
        p2036.get("result_kind")
        == "OPEN_NEW_PREMISE_BACKGROUND_TRANSPORT_CONTRACT_CANDIDATE_EXPORTED__NO_GLOBALIZATION"
        and as_bool(p2036_checks.get("candidate_contract_syntactically_complete"))
    )

    # First explicit same-scheme finite-part seed: identity correction at sigma2=0.
    delta_fp_same_scheme = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    map_same_scheme = [[I3[i][j] + delta_fp_same_scheme[i][j] for j in range(3)] for i in range(3)]
    residual_vs_identity = linf(mat_sub(map_same_scheme, I3))

    map_seed_exported = local_source_ready and obstruction_context_ready and contract_context_ready
    c2_basis_preserving_seed = map_seed_exported and residual_vs_identity == 0.0
    c3_theorem_proven = False

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2037",
        "stage_id": "S987",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_EXPLICIT_SAME_SCHEME_FINITE_PART_MAP_SEED_WITH_TRACE__C2_C3_STILL_OPEN"
            if map_seed_exported
            else "OPEN_SAME_SCHEME_FINITE_PART_MAP_SEED_BLOCKED"
        ),
        "route": "strict_only_with_local_quotient_transport_seed",
        "depends_on": {
            "p2034_present": p2034.get("_missing") is None,
            "p2035_present": p2035.get("_missing") is None,
            "p2036_present": p2036.get("_missing") is None,
            "p2034_local_source_ready": local_source_ready,
            "p2035_obstruction_context_ready": obstruction_context_ready,
            "p2036_contract_context_ready": contract_context_ready,
        },
        "input_hashes": {
            "p2034_json_sha256": file_sha256(GEN / "p2034_s984_strict_task1_quotient_only_renormalization_theorem.json"),
            "p2035_json_sha256": file_sha256(GEN / "p2035_s985_strict_task1_quotient_background_transport_obstruction_theorem.json"),
            "p2036_json_sha256": file_sha256(GEN / "p2036_s986_strict_task1_quotient_background_transport_candidate_contract.json"),
        },
        "same_scheme_finite_part_map_seed": {
            "map_id": "M_same_scheme_seed_v1",
            "basis": BASIS,
            "domain": "quotient coefficients (R2_bar,Ric2_bar,Riem2_bar)",
            "codomain": "quotient coefficients (R2_bar,Ric2_bar,Riem2_bar)",
            "sigma2_anchor": 0.0,
            "delta_finite_part_matrix": delta_fp_same_scheme,
            "transport_matrix": map_same_scheme,
            "residual_vs_identity_linf": residual_vs_identity,
            "scope": "same renormalization scheme seed only",
            "warning": "Seed map is not a theorem-grade finite-part transport proof.",
        },
        "c2_c3_gate_update": {
            "C2_basis_preserving_map_seed": "EXPORTED" if c2_basis_preserving_seed else "OPEN",
            "C3_finite_part_scheme_transport_theorem": "OPEN",
            "C2_C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_discharge": [
                "nontrivial anisotropic finite-part correction matrix from loop data",
                "proof of subtraction compatibility on same operator basis",
                "cross-background witness that transported targets agree",
            ],
        },
        "gatekeeper_checks": {
            "map_seed_exported": map_seed_exported,
            "basis_is_rank3_quotient_basis": BASIS == ["R2_bar", "Ric2_bar", "Riem2_bar"],
            "residual_vs_identity_zero": residual_vs_identity == 0.0,
            "c2_basis_preserving_seed_exported": c2_basis_preserving_seed,
            "c3_theorem_proven": c3_theorem_proven,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2037 S987: strict Task-1 same-scheme finite-part map seed",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported first explicit same-scheme finite-part seed map on quotient basis",
        "`R2_bar/Ric2_bar/Riem2_bar`: identity seed at `sigma2=0` with zero finite-part",
        "correction matrix.",
        "",
        "## Gate update",
        "",
        "- `C2`: seed exported (basis-preserving at anchor).",
        "- `C3`: theorem still open (no loop-derived correction + no subtraction-compatibility proof).",
        "",
        "## Discipline",
        "",
        "No tensor-component claims, no background-globalization claim, no ToE closure claim.",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
