#!/usr/bin/env python3
"""QW-2391: dual action-level anchor frontier alignment gate.

Compares anchor branch frontier (QW-2390) with existing foundational frontier
(QW-2296) to avoid false 'new closure' claims.
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
    q2390 = load("report_qw2390_dual_action_level_anchor_provider_execution_status_gate.json")
    q2296 = load("report_qw2296_dual_action_level_provider_execution_status_gate.json")

    curr_cut = norm_cut(q2390.get("blocker_cut", []))
    base_cut = norm_cut(q2296.get("blocker_cut", []))

    flags = {
        "q2390_execution_status_present": q2390.get("verdict")
        == "DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_FOUNDATIONAL_DERIVATION_SYMBOLS",
        "q2296_execution_status_present": q2296.get("verdict")
        == "DUAL_ACTION_LEVEL_PROVIDER_EXECUTION_STATUS_GATE_PASS_PARTIAL_BLOCKED_BY_MISSING_FOUNDATIONAL_DERIVATION_SYMBOLS",
        "current_blocker_cut_size_two": len(curr_cut) == 2,
        "baseline_blocker_cut_size_two": len(base_cut) == 2,
        "foundational_frontier_aligned": curr_cut == base_cut,
        "no_new_theorem_level_closure_claim": True,
        "execution_completed": False,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_ACTION_LEVEL_ANCHOR_FRONTIER_ALIGNMENT_GATE_PASS_ALIGNED_WITH_FOUNDATIONAL_CHAIN"
        if flags["q2390_execution_status_present"]
        and flags["q2296_execution_status_present"]
        and flags["current_blocker_cut_size_two"]
        and flags["baseline_blocker_cut_size_two"]
        and flags["foundational_frontier_aligned"]
        else "DUAL_ACTION_LEVEL_ANCHOR_FRONTIER_ALIGNMENT_GATE_FAIL"
    )

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2390": "report_qw2390_dual_action_level_anchor_provider_execution_status_gate.json",
            "q2296": "report_qw2296_dual_action_level_provider_execution_status_gate.json",
        },
        "current_blocker_cut_normalized": curr_cut,
        "baseline_blocker_cut_normalized": base_cut,
        "scope_boundary": {
            "frontier_alignment_established": curr_cut == base_cut,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2391_dual_action_level_anchor_frontier_alignment.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "current_blocker_cut": q2390.get("blocker_cut", []),
        "baseline_blocker_cut": q2296.get("blocker_cut", []),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "REUSE_FOUNDATIONAL_DERIVATION_CHAIN_FROM_QW2297_WITH_NONCYCLIC_GUARDS",
    }

    out_json = ROOT / "report_qw2391_dual_action_level_anchor_frontier_alignment_gate.json"
    out_md = ROOT / "RAPORT_QW2391_DUAL_ACTION_LEVEL_ANCHOR_FRONTIER_ALIGNMENT_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2391: DUAL ACTION LEVEL ANCHOR FRONTIER ALIGNMENT GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- foundational_frontier_aligned: `{flags['foundational_frontier_aligned']}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "aligned": flags["foundational_frontier_aligned"]}))


if __name__ == "__main__":
    main()
