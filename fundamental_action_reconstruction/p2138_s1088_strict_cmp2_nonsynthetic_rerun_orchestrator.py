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
OUT = GEN / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.json"
MD = GEN / "p2138_s1088_strict_cmp2_nonsynthetic_rerun_orchestrator.md"

SCHEMA_VERSION = "p2138_s1088_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"

PIPELINE = [
    ROOT / "p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py",
    ROOT / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.py",
    ROOT / "p2133_s1083_strict_cmp2_real_extension_merge_contract.py",
    ROOT / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.py",
    ROOT / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.py",
]


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def run_script(path: Path) -> tuple[bool, str]:
    try:
        subprocess.run([sys.executable, str(path)], check=True)
        return True, "OK"
    except subprocess.CalledProcessError as e:
        return False, f"FAILED({e.returncode})"


def main() -> None:
    GEN.mkdir(exist_ok=True)

    ext_payload = load(EXT) if EXT.exists() else {"rows": []}
    ext_rows = ext_payload.get("rows") or []

    execution_log: list[dict[str, Any]] = []
    for script in PIPELINE:
        ok, msg = run_script(script)
        execution_log.append({"script": script.name, "ok": ok, "message": msg})
        if not ok:
            break

    p2132 = load(GEN / "p2132_s1082_strict_cmp2_real_support_acquisition_gate.json")
    p2133 = load(GEN / "p2133_s1083_strict_cmp2_real_extension_merge_contract.json")
    p2134 = load(GEN / "p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json")
    p2135 = load(GEN / "p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json")

    k2132 = p2132.get("result_kind", "OPEN_MISSING")
    k2133 = p2133.get("result_kind", "OPEN_MISSING")
    k2134 = p2134.get("result_kind", "OPEN_MISSING")
    k2135 = p2135.get("result_kind", "OPEN_MISSING")

    fully_unblocked = (
        k2132 == "PASS_STRICT_CMP2_REAL_SUPPORT_ACQUISITION_GATE_WITH_TRACE"
        and k2133 == "PASS_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_WITH_TRACE"
        and k2134 == "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_COMPARISON_AUDIT_WITH_TRACE"
        and k2135 == "PASS_STRICT_CMP2_MERGED_REAL_BLOCK_VARIANT_STABILITY_AUDIT_WITH_TRACE"
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2138",
        "stage_id": "S1088",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_NONSYNTHETIC_RERUN_ORCHESTRATOR_WITH_TRACE" if fully_unblocked else "OPEN_STRICT_CMP2_NONSYNTHETIC_RERUN_ORCHESTRATOR_BLOCKED",
        "nonsynthetic_rerun_orchestrator": {
            "target_extension_file": str(EXT.relative_to(ROOT)),
            "rows_detected": len(ext_rows),
            "pipeline_execution_log": execution_log,
            "post_run_result_kinds": {
                "p2132": k2132,
                "p2133": k2133,
                "p2134": k2134,
                "p2135": k2135,
            },
            "fully_unblocked": fully_unblocked,
            "scope_limit": "orchestration/readout only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2139_candidate",
            "goal": "if fully unblocked, freeze non-synthetic CI inflation comparison report and move to next formal obstruction track",
        },
        "gatekeeper_checks": {
            "extension_rows_present": len(ext_rows) > 0,
            "p2132_pass": k2132.startswith("PASS_"),
            "p2133_pass": k2133.startswith("PASS_"),
            "p2134_pass": k2134.startswith("PASS_"),
            "p2135_pass": k2135.startswith("PASS_"),
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
            "# P2138 S1088: strict CMP2 non-synthetic rerun orchestrator",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Extension rows detected: `{len(ext_rows)}`",
            "",
            "This stage orchestrates P2132->P2135 rerun and reports whether non-synthetic path is fully unblocked.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
