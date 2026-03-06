#!/usr/bin/env python3
"""QW-2389: dual action-level anchor-provider discharge packet gate.

Builds execution packet for two core obligations from QW-2388.
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
    q2388 = load("report_qw2388_dual_action_level_anchor_provider_minimal_blocker_cut_gate.json")
    minimal_cut = q2388.get("minimal_blocker_cut", [])

    obligations = []
    for row in minimal_cut:
        obligations.append(
            {
                "id": row["id"],
                "branch": row["branch"],
                "provider_symbol": row["symbol"],
                "required_statement": row["required_statement"],
                "proof_kind": row["proof_kind"],
                "acceptance_criteria": [
                    "THEOREM_SYMBOL_EXISTS",
                    "NO_AXIOM_TOKENS_IN_PROVIDER_FILE",
                    "NO_DERIVED_OR_PENDING_TOKENS_IN_PROVIDER_FILE",
                    "LEAN_MACHINE_CHECK_EXECUTED",
                    "BLOCKER_CUT_IF_NONZERO_EXIT_IS_EXPLICIT",
                ],
            }
        )

    flags = {
        "q2388_minimal_cut_present": q2388.get("verdict")
        == "DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "minimal_cut_size_two": len(minimal_cut) == 2,
        "obligation_packet_size_two": len(obligations) == 2,
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    verdict = (
        "DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_DISCHARGE_PACKET_GATE_PASS_PACKET_READY"
        if flags["q2388_minimal_cut_present"] and flags["minimal_cut_size_two"] and flags["obligation_packet_size_two"]
        else "DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_DISCHARGE_PACKET_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2388_dual_action_level_anchor_provider_minimal_blocker_cut_gate.json",
        "obligations": obligations,
        "scope_boundary": {
            "packet_ready": True,
            "execution_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2389_dual_action_level_anchor_provider_discharge_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2388": "report_qw2388_dual_action_level_anchor_provider_minimal_blocker_cut_gate.json",
            "spec": spec_path.name,
        },
        "n_obligations": len(obligations),
    }
    proof_path = ROOT / "proof_object_qw2389_dual_action_level_anchor_provider_discharge_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "n_obligations": len(obligations),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EXECUTE_DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_DISCHARGE",
    }

    out_json = ROOT / "report_qw2389_dual_action_level_anchor_provider_discharge_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2389_DUAL_ACTION_LEVEL_ANCHOR_PROVIDER_DISCHARGE_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2389: DUAL ACTION LEVEL ANCHOR PROVIDER DISCHARGE PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_obligations: `{len(obligations)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations": len(obligations)}))


if __name__ == "__main__":
    main()
