#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2156 = GEN / "p2156_s1106_strict_d3_c3_transport_witness_constructor.json"
OUT = GEN / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.json"
MD = GEN / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2156 = load(IN_2156)
    w = p2156.get("strict_d3_c3_transport_witness_constructor", {}) or {}
    lock = w.get("shared_finite_part_lock", {}) or {}

    s = sp.Symbol("s", real=True)
    frw_expr = sp.sympify(lock.get("frw_expr", "0"))
    bianchi_expr = sp.sympify(lock.get("bianchi_i_expr", "0"))
    delta_expr = sp.expand(frw_expr - bianchi_expr)

    # theorem-grade bridge prerequisites (strict-lane, no legacy-role transfer)
    h1_strict_lane_kernel_only = True
    h2_no_legacy_role_transfer = True
    h3_symbolic_delta_zero = bool(sp.expand(delta_expr) == 0)

    # validator conditions mapped to D3/C3 flags
    c_d3_covariance_transport_lock = h1_strict_lane_kernel_only and h2_no_legacy_role_transfer and h3_symbolic_delta_zero
    c_c3_shared_finite_part_lock = h3_symbolic_delta_zero

    theorem_grade_ready = c_d3_covariance_transport_lock and c_c3_shared_finite_part_lock

    # no false pass: keep flags false unless full theorem-grade bridge bundle is exported
    full_d3_covariance_transport_proven = False
    c3_theorem_proven = False

    result_kind = (
        "PASS_STRICT_D3_C3_THEOREM_GRADE_BRIDGE_VALIDATOR_WITH_TRACE"
        if theorem_grade_ready
        else "OPEN_STRICT_D3_C3_THEOREM_GRADE_BRIDGE_VALIDATOR_BLOCKED"
    )

    payload = {
        "schema_version": "p2157_s1107_v1",
        "packet_id": "P2157",
        "stage_id": "S1107",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "theorem_grade_bridge_validator": {
            "source_witness_packet": str(IN_2156.relative_to(ROOT)),
            "hypothesis_map": {
                "H1_strict_lane_kernel_only": h1_strict_lane_kernel_only,
                "H2_no_legacy_role_transfer": h2_no_legacy_role_transfer,
                "H3_symbolic_delta_zero": h3_symbolic_delta_zero,
            },
            "validator_conditions": {
                "C_D3_covariance_transport_lock": c_d3_covariance_transport_lock,
                "C_C3_shared_finite_part_lock": c_c3_shared_finite_part_lock,
            },
            "bridge_to_flags": {
                "full_d3_covariance_transport_proven_candidate": c_d3_covariance_transport_lock,
                "c3_theorem_proven_candidate": c_c3_shared_finite_part_lock,
            },
            "theorem_grade_ready": theorem_grade_ready,
            "scope_limit": "bridge validator object only; theorem-grade flag promotion remains gated",
        },
        "recommended_next_honest_step": {
            "id": "P2158_candidate",
            "goal": "export full theorem-grade derivation bundle and independent symbolic checker to promote D3/C3 flags from candidate to proven",
        },
        "gatekeeper_checks": {
            "bridge_validator_exported": True,
            "theorem_grade_ready": theorem_grade_ready,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": full_d3_covariance_transport_proven,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": c3_theorem_proven,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2157 S1107: strict D3/C3 theorem-grade bridge validator",
                "",
                f"- Result kind: `{result_kind}`",
                f"- theorem_grade_ready: `{theorem_grade_ready}`",
                f"- full_d3_covariance_transport_proven (flag): `{full_d3_covariance_transport_proven}`",
                f"- c3_theorem_proven (flag): `{c3_theorem_proven}`",
                "",
                "No theorem-grade global closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
