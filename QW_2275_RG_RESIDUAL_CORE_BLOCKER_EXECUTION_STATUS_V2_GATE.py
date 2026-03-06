#!/usr/bin/env python3
"""QW-2275: RG residual core-blocker execution status v2 gate.

Upgrades execution status with strict non-axiomatic evidence filter.
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
    q2271 = load("report_qw2271_rg_residual_core_blocker_execution_status_gate.json")
    q2273 = load("report_qw2273_rg_residual_non_axiomatic_provider_evidence_gate.json")
    p2269 = load("proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json")

    n_total = int(q2271.get("n_obligations_total", 0))
    strict_candidates = int(q2273.get("n_strict_non_axiomatic_candidates", 0))
    n_satisfied = 1 if strict_candidates > 0 and n_total == 1 else 0

    flags = {
        "q2271_execution_status_present": q2271.get("verdict")
        == "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_GATE_PASS_PARTIAL_PENDING",
        "q2273_evidence_present": q2273.get("verdict")
        == "RG_RESIDUAL_NON_AXIOMATIC_PROVIDER_EVIDENCE_GATE_PASS_PARTIAL_NO_STRICT_CANDIDATE",
        "single_residual_obligation_present": int(p2269.get("n_obligations", 0)) == 1,
        "strict_non_axiomatic_provider_found": strict_candidates > 0,
        "all_obligations_satisfied_strict": n_satisfied == n_total and n_total > 0,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V2_GATE_PASS_PARTIAL_PENDING_STRICT_NON_AXIOMATIC"
        if flags["q2271_execution_status_present"] and flags["q2273_evidence_present"]
        else "RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V2_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2271": "report_qw2271_rg_residual_core_blocker_execution_status_gate.json",
            "q2273": "report_qw2273_rg_residual_non_axiomatic_provider_evidence_gate.json",
            "p2269": "proof_object_qw2269_rg_residual_core_blocker_discharge_spec.json",
        },
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict": n_satisfied,
        "n_strict_non_axiomatic_candidates": strict_candidates,
    }

    proof_path = ROOT / "proof_object_qw2275_rg_residual_core_blocker_execution_status_v2.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations_total": n_total,
        "n_obligations_satisfied_strict": n_satisfied,
        "n_strict_non_axiomatic_candidates": strict_candidates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "PROVE_RG_RESIDUAL_O1_WITH_STRICT_NON_AXIOMATIC_PROVIDER",
    }

    out_json = ROOT / "report_qw2275_rg_residual_core_blocker_execution_status_v2_gate.json"
    out_md = ROOT / "RAPORT_QW2275_RG_RESIDUAL_CORE_BLOCKER_EXECUTION_STATUS_V2_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2275: RG RESIDUAL CORE BLOCKER EXECUTION STATUS V2 GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- obligations satisfied (strict): `{n_satisfied}/{n_total}`",
                f"- n_strict_non_axiomatic_candidates: `{strict_candidates}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations_satisfied_strict": n_satisfied, "n_obligations_total": n_total}))


if __name__ == "__main__":
    main()
