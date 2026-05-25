#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2094 = GEN / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.json"
IN_2095 = GEN / "p2095_s1045_strict_b1_gb_derived_channel_certificate.json"
OUT = GEN / "p2096_s1046_strict_b1_renormalization_closure_contract.json"
MD = GEN / "p2096_s1046_strict_b1_renormalization_closure_contract.md"

SCHEMA_VERSION = "p2096_s1046_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2094 = load(IN_2094)
    p2095 = load(IN_2095)

    checks_2094 = p2094.get("gatekeeper_checks", {})
    checks_2095 = p2095.get("gatekeeper_checks", {})

    pre_ready = (
        checks_2094.get("quotient_3x3_residual_small") is True
        and checks_2094.get("gb_row_consistent_as_derived_channel") is True
        and checks_2095.get("gb_derived_identity_operator_level") is True
    )

    q = p2094.get("quotient_rank_repair", {})
    sol = q.get("solution", {})
    closure_rows = [
        {
            "delta_label": "delta_c_gr_1",
            "channel": "R2",
            "coefficient": sol.get("a_R2"),
        },
        {
            "delta_label": "delta_c_gr_2",
            "channel": "Ric2",
            "coefficient": sol.get("a_Ric2"),
        },
        {
            "delta_label": "delta_c_gr_3",
            "channel": "Riem2",
            "coefficient": sol.get("a_Riem2"),
        },
        {
            "delta_label": "delta_c_gr_4",
            "channel": "GB (derived)",
            "coefficient": "derived via GB=Riem2-4Ric2+R2; no independent solve",
        },
    ]

    contract_conditions = {
        "C_RANK": "3x3 quotient solve residual <= tolerance",
        "C_GB_DERIVED": "operator-level and RHS-level GB derived identity residual <= tolerance",
        "C_SCOPE": "contract valid only for current strict B1 projection basis and exported conventions",
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2096",
        "stage_id": "S1046",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_STRICT_B1_RENORMALIZATION_CLOSURE_CONTRACT_WITH_TRACE__QUOTIENT_SCOPE_ONLY"
            if pre_ready
            else "OPEN_STRICT_B1_RENORMALIZATION_CLOSURE_CONTRACT_BLOCKED"
        ),
        "depends_on": {
            "p2094_present": p2094.get("_missing") is None,
            "p2095_present": p2095.get("_missing") is None,
            "preconditions_ready": pre_ready,
        },
        "closure_contract": {
            "counterterm_rows": closure_rows,
            "contract_conditions": contract_conditions,
            "residuals": {
                "quotient_residual_abs_max": q.get("residual_abs_max"),
                "gb_derived_operator_residual_abs_max": (p2095.get("derived_channel_certificate", {}) or {}).get("operator_residual_abs_max"),
                "gb_derived_rhs_residual": (p2095.get("derived_channel_certificate", {}) or {}).get("rhs_residual"),
            },
            "scope_limit": "Not a full 4-channel independence theorem; not global background-family closure.",
        },
        "c3_gate_update": {
            "C3_strict_b1_renormalization_closure_contract": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": pre_ready,
            "quotient_scope_contract_exported": True,
            "full_4channel_independence_proven": False,
            "global_background_family_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2096 S1046: strict B1 renormalization closure contract",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Quotient residual max: `{payload['closure_contract']['residuals']['quotient_residual_abs_max']}`",
            f"- GB derived operator residual max: `{payload['closure_contract']['residuals']['gb_derived_operator_residual_abs_max']}`",
            "",
            "This stage exports a strict B1 renormalization closure contract in quotient scope only.",
            "No full 4-channel independence or global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
