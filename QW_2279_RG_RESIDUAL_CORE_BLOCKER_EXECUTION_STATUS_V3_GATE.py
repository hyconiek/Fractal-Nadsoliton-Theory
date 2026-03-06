#!/usr/bin/env python3
"""QW-2279: RG residual core-blocker execution status v3 gate.

Combines lexical strict-candidate status (QW-2275) with machine-checkable
construction attempt status (QW-2277).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> None:
    q2275 = load("report_qw2275_rg_residual_core_blocker_execution_status_v2_gate.json")
    q2277 = load("report_qw2277_rg_residual_strict_non_axiomatic_provider_construction_gate.json")
    p2269 = load("proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json")

    n_total = int(q2275.get("n_obligations_total", 0))
    n_lexical = int(q2275.get("n_strict_non_axiomatic_candidates", 0))
    machine_exit = q2277.get("machine_check_exit_code")
    machine_success = machine_exit == 0

    n_satisfied_strict = 1 if (n_total == 1 and n_lexical > 0 and machine_success) else 0

    flags = {
        "q2275_v2_present": q2275.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_PENDING_STRICT_NON_AXIOMATIC",
        "q2277_machine_attempt_present": q2277.get("verdict")
        in {
            "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_PASS_PARTIAL_OBSTRUCTION_CONFIRMED",
            "RG_RESIDUAL_STRICT_NON_AXIOMATIC_PROVIDER_CONSTRUCTION_GATE_PASS_STRICT_PROVIDER_CONSTRUCTED",
        },
        "single_residual_obligation_present": int(p2269.get("n_obligations", 0)) == 1,
        "lexical_strict_candidate_found": n_lexical > 0,
        "machine_checkable_provider_constructed": machine_success,
        "all_obligations_satisfied_strict_v3": n_satisfied_strict == n_total and n_total > 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V3_GATE_PASS_PARTIAL_PENDING_MACHINE_CHECKABLE_NON_AXIOMATIC"
        if flags["q2275_v2_present"] and flags["q2277_machine_attempt_present"]
        else "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V3_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2275": "report_qw2275_rg_residual_core_blocker_execution_status_v2_gate.json",
            "q2277": "report_qw2277_rg_residual_strict_non_axiomatic_provider_construction_gate.json",
            "p2269": "proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json",
        },
        "n_obligations_total": n_total,
        "n_strict_non_axiomatic_candidates_lexical": n_lexical,
        "machine_check_exit_code": machine_exit,
        "n_obligations_satisfied_strict_v3": n_satisfied_strict,
    }

    proof_path = ROOT / "proof_object_qw2279_rg_residual_core_blocker_execution_status_v3.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations_total": n_total,
        "n_strict_non_axiomatic_candidates_lexical": n_lexical,
        "machine_check_exit_code": machine_exit,
        "n_obligations_satisfied_strict_v3": n_satisfied_strict,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_RG_PROVIDER_TO_MACHINE_CHECKABLE_NON_AXIOMATIC_THEOREM"
            if n_satisfied_strict < n_total
            else "UPDATE_RG_SCOPE_CLOSURE_STATUS"
        ),
    }

    out_json = ROOT / "report_qw2279_rg_residual_core_blocker_execution_status_v3_gate.json"
    out_md = ROOT / "RAPORT_QW2279_RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V3_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2279: RG RESIDUAL CORE BLOCKER EXECUTION STATUS V3 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- strict lexical candidates: `{n_lexical}`",
                f"- machine_check_exit_code: `{machine_exit}`",
                f"- obligations satisfied (strict v3): `{n_satisfied_strict}/{n_total}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "verdict": verdict,
                "n_obligations_satisfied_strict_v3": n_satisfied_strict,
                "n_obligations_total": n_total,
            }
        )
    )


if __name__ == "__main__":
    main()
