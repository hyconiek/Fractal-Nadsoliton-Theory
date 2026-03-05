#!/usr/bin/env python3
"""QW-2250: QFT export axiom-free obligation packet gate.

Builds machine-checkable obligation packet for missing non-axiomatic QFT export theorem.
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
    q2244 = load("report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json")
    q2246 = load("report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json")
    q2248 = load("report_qw2248_qft_export_axiomatic_dependency_gate.json")

    obligations = {
        "QFT_EXPORT_O1": {
            "statement": "theorem QFT_CanonicalAction_to_Positivity_EXPORT : (FINActionComplete ∧ ConstructiveNonPerturbativeScheme) -> PositivityToReconstruction",
            "required": True,
            "acceptance": "symbol exists in *.lean",
        },
        "QFT_EXPORT_O2": {
            "statement": "proof of QFT_CanonicalAction_to_Positivity_EXPORT uses no *DerivedOrPending symbol",
            "required": True,
            "acceptance": "static proof refs contain zero tokens ending with _DerivedOrPending",
        },
        "QFT_EXPORT_O3": {
            "statement": "proof of QFT_CanonicalAction_to_Positivity_EXPORT is axiom-free",
            "required": True,
            "acceptance": "no local/global axiom dependency in theorem dependency chain",
        },
        "QFT_EXPORT_O4": {
            "statement": "QW-2244 rerun passes with canonical_export_symbol_exists=True",
            "required": True,
            "acceptance": "report_qw2244... flags.canonical_export_symbol_exists == true",
        },
    }

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2244": "report_qw2244_qft_dax1_non_axiomatic_provider_attempt_gate.json",
            "q2246": "report_qw2246_qft_dax1_axiom_free_candidate_scan_gate.json",
            "q2248": "report_qw2248_qft_export_axiomatic_dependency_gate.json",
        },
        "obligations": obligations,
        "current_state": {
            "q2244_verdict": q2244.get("verdict"),
            "q2246_verdict": q2246.get("verdict"),
            "q2248_verdict": q2248.get("verdict"),
            "canonical_export_symbol_exists": q2244.get("flags", {}).get("canonical_export_symbol_exists", False),
            "axiom_free_candidates_found": q2246.get("scan_summary", {}).get("n_axiom_free_candidates", 0),
            "axiom_dependency_certified": q2248.get("flags", {}).get("dependency_chain_hits_axiom_layer", False),
        },
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2250_qft_export_axiom_free_obligation_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2244_attempt_blocker_confirmed": q2244.get("verdict")
        == "QFT_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "q2246_candidate_scan_confirmed": q2246.get("verdict")
        == "QFT_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_PASS_PARTIAL_NO_AXIOM_FREE_CANDIDATE",
        "q2248_dependency_certified": q2248.get("verdict")
        == "QFT_EXPORT_AXIOMATIC_DEPENDENCY_GATE_PASS_PARTIAL_AXIOM_FREE_EXPORT_ABSENT",
        "obligation_packet_written": proof_path.exists(),
        "all_obligations_explicit": len(obligations) == 4,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2244_attempt_blocker_confirmed"]
        and flags["q2246_candidate_scan_confirmed"]
        and flags["q2248_dependency_certified"]
        and flags["obligation_packet_written"]
        and flags["all_obligations_explicit"]
    )

    verdict = (
        "QFT_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE_PASS_PACKET_READY_EXPORT_PENDING"
        if core_ok
        else "QFT_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "obligation_ids": list(obligations.keys()),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "CONSTRUCT_NON_AXIOMATIC_QFT_EXPORT_THEOREM_AND_RERUN_QW2244_QW2246_QW2248",
    }

    out_json = ROOT / "report_qw2250_qft_export_axiom_free_obligation_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2250_QFT_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2250: QFT EXPORT AXIOM-FREE OBLIGATION PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zbudowano jawny pakiet 4 obligacji formalnych dla export theorem QFT.",
                "- Pakiet gotowy do wykonania; export theorem nadal pending.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2250_qft_export_axiom_free_obligation_packet_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
