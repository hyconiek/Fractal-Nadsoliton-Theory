#!/usr/bin/env python3
"""QW-2403: dual kernel-operator anchor frontier alignment gate.

Compares kernel-operator anchor frontier (QW-2402) with historical
kernel-operator execution frontier (QW-2305) to avoid false 'new closure'
claims.
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


def norm_cut(rows: list[dict[str, Any]]) -> list[tuple[str, str]]:
    return sorted((str(r.get("branch", "")), str(r.get("symbol", ""))) for r in rows)


def main() -> None:
    q2402 = load("report_qw2402_dual_kernel_operator_anchor_execution_status_gate.json")
    q2305 = load("report_qw2305_dual_kernel_operator_closure_execution_status_gate.json")

    curr = norm_cut(q2402.get("blocker_cut", []))
    base = norm_cut(q2305.get("blocker_cut", []))

    flags = {
        "q2402_execution_status_present": q2402.get("verdict")
        == "DUAL_KERNEL_OPERATOR_ANCHOR_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_SPECTRAL_CLOSURE_THEOREMS",
        "q2305_execution_status_present": q2305.get("verdict")
        == "DUAL_KERNEL_OPERATOR_CLOSURE_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_KERNEL_SPECTRAL_CLOSURE_THEOREMS",
        "current_blocker_cut_size_two": len(curr) == 2,
        "baseline_blocker_cut_size_two": len(base) == 2,
        "kernel_spectral_frontier_aligned": curr == base,
        "no_new_theorem_level_closure_claim": True,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_KERNEL_OPERATOR_ANCHOR_FRONTIER_ALIGNMENT_GATE_PASS_ALIGNED_WITH_KERNEL_SPECTRAL_CHAIN"
        if flags["q2402_execution_status_present"]
        and flags["q2305_execution_status_present"]
        and flags["current_blocker_cut_size_two"]
        and flags["baseline_blocker_cut_size_two"]
        and flags["kernel_spectral_frontier_aligned"]
        else "DUAL_KERNEL_OPERATOR_ANCHOR_FRONTIER_ALIGNMENT_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2402": "report_qw2402_dual_kernel_operator_anchor_execution_status_gate.json",
            "q2305": "report_qw2305_dual_kernel_operator_closure_execution_status_gate.json",
        },
        "current_blocker_cut_normalized": curr,
        "baseline_blocker_cut_normalized": base,
        "scope_boundary": {
            "frontier_alignment_established": curr == base,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2403_dual_kernel_operator_anchor_frontier_alignment.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "current_blocker_cut": q2402.get("blocker_cut", []),
        "baseline_blocker_cut": q2305.get("blocker_cut", []),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "REUSE_KERNEL_SPECTRAL_CLOSURE_CHAIN_FROM_QW2306_WITH_NONCYCLIC_GUARDS",
    }

    out_json = ROOT / "report_qw2403_dual_kernel_operator_anchor_frontier_alignment_gate.json"
    out_md = ROOT / "RAPORT_QW2403_DUAL_KERNEL_OPERATOR_ANCHOR_FRONTIER_ALIGNMENT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2403: DUAL KERNEL OPERATOR ANCHOR FRONTIER ALIGNMENT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- kernel_spectral_frontier_aligned: `{flags['kernel_spectral_frontier_aligned']}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "aligned": flags["kernel_spectral_frontier_aligned"]}))


if __name__ == "__main__":
    main()
