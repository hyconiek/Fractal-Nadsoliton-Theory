#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
EXT = GEN / "cmp2_backend_rows_extension_v1.json"
TEMPLATE = GEN / "cmp2_backend_rows_extension_v1.template.json"
OUT = GEN / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json"
MD = GEN / "p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.md"

SCHEMA_VERSION = "p2142_s1092_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"

PIPELINE = [
    ROOT / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.py",
    ROOT / "p2133_s1083_strict_cmp2_real_extension_merge_contract.py",
    ROOT / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.py",
    ROOT / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.py",
    ROOT / "p2141_s1091_strict_cmp2_flag_outcome_matrix.py",
]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def run(script: Path) -> tuple[bool, str]:
    try:
        subprocess.run([sys.executable, str(script)], check=True)
        return True, "OK"
    except subprocess.CalledProcessError as e:
        return False, f"FAILED({e.returncode})"


def main() -> None:
    GEN.mkdir(exist_ok=True)

    ext = load(EXT) if EXT.exists() else {"rows": []}
    rows = ext.get("rows") or []

    log = []
    for s in PIPELINE:
        ok, msg = run(s)
        log.append({"script": s.name, "ok": ok, "message": msg})
        if not ok:
            break

    p2132 = load(GEN / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.json")
    p2133 = load(GEN / "p2133_s1083_strict_cmp2_real_extension_merge_contract.json")
    p2134 = load(GEN / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json")
    p2135 = load(GEN / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json")
    p2141 = load(GEN / "p2141_s1091_strict_cmp2_flag_outcome_matrix.json")

    rk = {
        "p2132": p2132.get("result_kind", "OPEN_MISSING"),
        "p2133": p2133.get("result_kind", "OPEN_MISSING"),
        "p2134": p2134.get("result_kind", "OPEN_MISSING"),
        "p2135": p2135.get("result_kind", "OPEN_MISSING"),
        "p2141": p2141.get("result_kind", "OPEN_MISSING"),
    }

    fully_ready = all(rk[k].startswith("PASS_") for k in ["p2132", "p2133", "p2134", "p2135"])

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2142",
        "stage_id": "S1092",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_REAL_DATA_HANDOFF_AND_RERUN_EXECUTOR_WITH_TRACE" if fully_ready else "OPEN_STRICT_CMP2_REAL_DATA_HANDOFF_AND_RERUN_EXECUTOR_BLOCKED",
        "real_data_handoff_and_rerun_executor": {
            "extension_file": str(EXT.relative_to(ROOT)),
            "template_file": str(TEMPLATE.relative_to(ROOT)),
            "rows_detected": len(rows),
            "pipeline_log": log,
            "post_rerun_result_kinds": rk,
            "fully_ready_for_real_ci_stability_assessment": fully_ready,
            "scope_limit": "execution/handoff status only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2143_candidate",
            "goal": "if blocked: deliver real extension rows and rerun P2142; if ready: perform final non-synthetic CI stability interpretation memo",
        },
        "gatekeeper_checks": {
            "extension_rows_present": len(rows) > 0,
            "p2132_pass": rk["p2132"].startswith("PASS_"),
            "p2133_pass": rk["p2133"].startswith("PASS_"),
            "p2134_pass": rk["p2134"].startswith("PASS_"),
            "p2135_pass": rk["p2135"].startswith("PASS_"),
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2142 S1092: strict CMP2 real-data handoff and rerun executor",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Rows detected: `{len(rows)}`",
            "",
            "This stage executes the requested P2132/P2133/P2134/P2135 + P2141 sequence and reports readiness for real CI stability assessment.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
