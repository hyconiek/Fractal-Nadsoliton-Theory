#!/usr/bin/env python3
"""QW-2292: dual physical-premise discharge packet gate.

Builds an explicit discharge packet for the two remaining physical bridge
premises identified by QW-2291.
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
    q2291 = load("report_qw2291_dual_single_premise_frontier_gate.json")
    frontier = q2291.get("remaining_frontier", [])

    obligations = []
    for row in frontier:
        branch = row.get("branch")
        premise = row.get("premise_symbol")
        stmt = row.get("statement")
        obligations.append(
            {
                "id": f"{branch}_PHYSICAL_PREMISE_O1",
                "branch": branch,
                "premise_symbol": premise,
                "target_statement": stmt,
                "required_outcome": {
                    "construct_theorem_symbol": premise,
                    "proof_kind": "nonlogical_physical_derivation",
                    "scope": "axiom-token-free machine-checkable theorem object",
                },
                "acceptance_criteria": [
                    "THEOREM_SYMBOL_EXISTS",
                    "NO_AXIOM_TOKENS_IN_PROVIDER_FILE",
                    "NO_DERIVED_OR_PENDING_TOKENS_IN_PROVIDER_FILE",
                    "LEAN_MACHINE_CHECK_EXIT_ZERO",
                    "ACTION_LEVEL_PREMISE_TRACEABILITY_EXPLICIT",
                ],
            }
        )

    flags = {
        "q2291_frontier_present": q2291.get("verdict")
        == "DUAL_SINGLE_PREMISE_FRONTIER_GATE_PASS_PARTIAL_FRONTIER_EXPLICIT",
        "frontier_size_two": len(frontier) == 2,
        "obligation_packet_size_two": len(obligations) == 2,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_PHYSICAL_PREMISE_DISCHARGE_PACKET_GATE_PASS_PACKET_READY"
        if flags["q2291_frontier_present"] and flags["frontier_size_two"] and flags["obligation_packet_size_two"]
        else "DUAL_PHYSICAL_PREMISE_DISCHARGE_PACKET_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2291_dual_single_premise_frontier_gate.json",
        "remaining_frontier": frontier,
        "obligations": obligations,
        "scope_boundary": {
            "packet_ready": True,
            "nonlogical_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2292_dual_physical_premise_discharge_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {"q2291": "report_qw2291_dual_single_premise_frontier_gate.json", "spec": spec_path.name},
        "n_frontier_items": len(frontier),
        "n_obligations": len(obligations),
    }
    proof_path = ROOT / "proof_object_qw2292_dual_physical_premise_discharge_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_frontier_items": len(frontier),
        "n_obligations": len(obligations),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EXECUTE_DUAL_PHYSICAL_PREMISE_NONLOGICAL_DISCHARGE",
    }

    out_json = ROOT / "report_qw2292_dual_physical_premise_discharge_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2292_DUAL_PHYSICAL_PREMISE_DISCHARGE_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2292: DUAL PHYSICAL PREMISE DISCHARGE PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_frontier_items: `{len(frontier)}`",
                f"- n_obligations: `{len(obligations)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations": len(obligations)}))


if __name__ == "__main__":
    main()
