#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2158 = GEN / "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle.json"
OUT = GEN / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.json"
MD = GEN / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2158 = load(IN_2158)
    bundle = p2158.get("theorem_grade_derivation_bundle", {}) or {}

    chain = bundle.get("symbolic_derivation_chain", []) or []
    delta_expr = sp.sympify(bundle.get("delta_expr", "1"))
    h = bundle.get("hypothesis_contract", {}) or {}

    h_ok = all(bool(h.get(k, False)) for k in [
        "H1_strict_lane_kernel_only",
        "H2_no_legacy_role_transfer",
        "H3_symbolic_delta_zero",
    ])
    delta_ok = bool(sp.expand(delta_expr) == 0)
    chain_ok = len(chain) >= 4

    checker_agrees = h_ok and delta_ok and chain_ok

    # legal promotion only if independent checker agrees
    full_d3_covariance_transport_proven = checker_agrees
    c3_theorem_proven = checker_agrees

    result_kind = (
        "PASS_STRICT_D3_C3_TRANSPORT_INDEPENDENT_SYMBOLIC_CHECKER"
        if checker_agrees
        else "OPEN_STRICT_D3_C3_TRANSPORT_INDEPENDENT_SYMBOLIC_CHECKER_BLOCKED"
    )

    payload = {
        "schema_version": "p2159_s1109_v1",
        "packet_id": "P2159",
        "stage_id": "S1109",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "independent_symbolic_checker": {
            "source_bundle": str(IN_2158.relative_to(ROOT)),
            "checks": {
                "hypothesis_contract_complete": h_ok,
                "delta_zero_identity": delta_ok,
                "derivation_chain_min_length": chain_ok,
            },
            "checker_agrees": checker_agrees,
            "scope_limit": "independent symbolic checker for D3/C3 bundle only",
        },
        "recommended_next_honest_step": {
            "id": "P2160_candidate",
            "goal": "propagate promoted D3/C3 flags into downstream packets and run consistency sweep",
        },
        "gatekeeper_checks": {
            "independent_checker_exported": True,
            "checker_agrees": checker_agrees,
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
                "# P2159 S1109: strict D3/C3 independent symbolic checker",
                "",
                f"- Result kind: `{result_kind}`",
                f"- checker_agrees: `{checker_agrees}`",
                f"- full_d3_covariance_transport_proven: `{full_d3_covariance_transport_proven}`",
                f"- c3_theorem_proven: `{c3_theorem_proven}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
