#!/usr/bin/env python3
"""QW-2249: RG export axiom-free obligation packet gate.

Builds machine-checkable obligation packet for missing non-axiomatic RG export theorem.
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
    q2243 = load("report_qw2243_rg_dax1_non_axiomatic_provider_attempt_gate.json")
    q2245 = load("report_qw2245_rg_dax1_axiom_free_candidate_scan_gate.json")
    q2247 = load("report_qw2247_rg_export_axiomatic_dependency_gate.json")

    obligations = {
        "RG_EXPORT_O1": {
            "statement": "theorem RG_CanonicalAction_to_WellPosedness_EXPORT : (FINActionComplete ∧ RGConstructiveMap) -> RGGlobalWellPosednessAllScales",
            "required": True,
            "acceptance": "symbol exists in *.lean",
        },
        "RG_EXPORT_O2": {
            "statement": "proof of RG_CanonicalAction_to_WellPosedness_EXPORT uses no *DerivedOrPending symbol",
            "required": True,
            "acceptance": "static proof refs contain zero tokens ending with _DerivedOrPending",
        },
        "RG_EXPORT_O3": {
            "statement": "proof of RG_CanonicalAction_to_WellPosedness_EXPORT is axiom-free",
            "required": True,
            "acceptance": "no local/global axiom dependency in theorem dependency chain",
        },
        "RG_EXPORT_O4": {
            "statement": "QW-2243 rerun passes with canonical_export_symbol_exists=True",
            "required": True,
            "acceptance": "report_qw2243... flags.canonical_export_symbol_exists == true",
        },
    }

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2243": "report_qw2243_rg_dax1_non_axiomatic_provider_attempt_gate.json",
            "q2245": "report_qw2245_rg_dax1_axiom_free_candidate_scan_gate.json",
            "q2247": "report_qw2247_rg_export_axiomatic_dependency_gate.json",
        },
        "obligations": obligations,
        "current_state": {
            "q2243_verdict": q2243.get("verdict"),
            "q2245_verdict": q2245.get("verdict"),
            "q2247_verdict": q2247.get("verdict"),
            "canonical_export_symbol_exists": q2243.get("flags", {}).get("canonical_export_symbol_exists", False),
            "axiom_free_candidates_found": q2245.get("scan_summary", {}).get("n_axiom_free_candidates", 0),
            "axiom_dependency_certified": q2247.get("flags", {}).get("dependency_chain_hits_axiom_layer", False),
        },
        "scope_boundary": {
            "dax1_non_axiomatic_provider_completed": False,
            "c1_theorem_discharge_completed": False,
            "o1c_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }

    proof_path = ROOT / "proof_object_qw2249_rg_export_axiom_free_obligation_packet.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    flags = {
        "q2243_attempt_blocker_confirmed": q2243.get("verdict")
        == "RG_DAX1_NON_AXIOMATIC_PROVIDER_ATTEMPT_GATE_PASS_PARTIAL_CANONICAL_EXPORT_MISSING",
        "q2245_candidate_scan_confirmed": q2245.get("verdict")
        == "RG_DAX1_AXIOM_FREE_CANDIDATE_SCAN_GATE_PASS_PARTIAL_NO_AXIOM_FREE_CANDIDATE",
        "q2247_dependency_certified": q2247.get("verdict")
        == "RG_EXPORT_AXIOMATIC_DEPENDENCY_GATE_PASS_PARTIAL_AXIOM_FREE_EXPORT_ABSENT",
        "obligation_packet_written": proof_path.exists(),
        "all_obligations_explicit": len(obligations) == 4,
        "dax1_non_axiomatic_provider_completed": False,
        "o1c_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2243_attempt_blocker_confirmed"]
        and flags["q2245_candidate_scan_confirmed"]
        and flags["q2247_dependency_certified"]
        and flags["obligation_packet_written"]
        and flags["all_obligations_explicit"]
    )

    verdict = (
        "RG_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE_PASS_PACKET_READY_EXPORT_PENDING"
        if core_ok
        else "RG_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE_FAIL"
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
        "required_next_step": "CONSTRUCT_NON_AXIOMATIC_RG_EXPORT_THEOREM_AND_RERUN_QW2243_QW2245_QW2247",
    }

    out_json = ROOT / "report_qw2249_rg_export_axiom_free_obligation_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2249_RG_EXPORT_AXIOM_FREE_OBLIGATION_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2249: RG EXPORT AXIOM-FREE OBLIGATION PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zbudowano jawny pakiet 4 obligacji formalnych dla export theorem RG.",
                "- Pakiet gotowy do wykonania; export theorem nadal pending.",
                "",
                "## Artifacts",
                "- JSON: `report_qw2249_rg_export_axiom_free_obligation_packet_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
