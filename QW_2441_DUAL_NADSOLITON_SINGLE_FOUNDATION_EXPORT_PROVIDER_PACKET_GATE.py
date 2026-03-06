#!/usr/bin/env python3
"""QW-2441: dual Nadsoliton single-foundation export-provider packet gate.

Builds a minimal two-obligation packet under the ontological assumption:
Nadsoliton is the single foundation, RG/QFT layers are emergent.
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
    q2440 = load("report_qw2440_grep_frontier_single_foundation_audit_gate.json")

    obligations = [
        {
            "id": "RG_NADSOLITON_SINGLE_FOUNDATION_O1",
            "branch": "L12",
            "foundation_postulate": "NadsolitonSingleFoundation",
            "target_export_symbol": "RG_CanonicalAction_to_WellPosedness_EXPORT",
            "required_provider_symbol": "RG_NadsolitonSingleFoundationToWellPosedness_EXPORT",
            "target_statement": "theorem RG_NadsolitonSingleFoundationToWellPosedness_EXPORT : NadsolitonSingleFoundation -> RG_CanonicalAction_to_WellPosedness_EXPORT",
            "acceptance_criteria": [
                "THEOREM_SYMBOL_EXISTS",
                "NO_AXIOM_TOKENS_IN_PROVIDER_FILE",
                "NO_DERIVED_OR_PENDING_TOKENS_IN_PROVIDER_FILE",
                "LEAN_MACHINE_CHECK_EXIT_ZERO",
                "FOUNDATION_TO_EXPORT_TRACEABILITY_EXPLICIT",
            ],
        },
        {
            "id": "QFT_NADSOLITON_SINGLE_FOUNDATION_O1",
            "branch": "L5",
            "foundation_postulate": "NadsolitonSingleFoundation",
            "target_export_symbol": "QFT_CanonicalAction_to_Positivity_EXPORT",
            "required_provider_symbol": "QFT_NadsolitonSingleFoundationToPositivity_EXPORT",
            "target_statement": "theorem QFT_NadsolitonSingleFoundationToPositivity_EXPORT : NadsolitonSingleFoundation -> QFT_CanonicalAction_to_Positivity_EXPORT",
            "acceptance_criteria": [
                "THEOREM_SYMBOL_EXISTS",
                "NO_AXIOM_TOKENS_IN_PROVIDER_FILE",
                "NO_DERIVED_OR_PENDING_TOKENS_IN_PROVIDER_FILE",
                "LEAN_MACHINE_CHECK_EXIT_ZERO",
                "FOUNDATION_TO_EXPORT_TRACEABILITY_EXPLICIT",
            ],
        },
    ]

    flags = {
        "q2440_audit_pass_with_blockers_explicit": q2440.get("verdict")
        == "GREP_FRONTIER_SINGLE_FOUNDATION_AUDIT_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
        "dual_packet_size_two": len(obligations) == 2,
        "single_foundation_assumption_explicit": all(o["foundation_postulate"] == "NadsolitonSingleFoundation" for o in obligations),
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    if flags["q2440_audit_pass_with_blockers_explicit"] and flags["dual_packet_size_two"]:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_EXPORT_PROVIDER_PACKET_GATE_PASS_PACKET_READY"
        required_next_step = "EXECUTE_DUAL_NADSOLITON_SINGLE_FOUNDATION_PROVIDER_ATTEMPT"
    else:
        verdict = "DUAL_NADSOLITON_SINGLE_FOUNDATION_EXPORT_PROVIDER_PACKET_GATE_FAIL"
        required_next_step = "REPAIR_SINGLE_FOUNDATION_PACKET_INPUTS"

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2440_grep_frontier_single_foundation_audit_gate.json",
        "obligations": obligations,
        "scope_boundary": {
            "packet_ready": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2441_dual_nadsoliton_single_foundation_export_provider_packet.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2440": "report_qw2440_grep_frontier_single_foundation_audit_gate.json",
            "spec": spec_path.name,
        },
        "n_obligations": len(obligations),
    }
    proof_path = ROOT / "proof_object_qw2441_dual_nadsoliton_single_foundation_export_provider_packet.json"
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
        "required_next_step": required_next_step,
    }

    out_json = ROOT / "report_qw2441_dual_nadsoliton_single_foundation_export_provider_packet_gate.json"
    out_md = ROOT / "RAPORT_QW2441_DUAL_NADSOLITON_SINGLE_FOUNDATION_EXPORT_PROVIDER_PACKET_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2441: DUAL NADSOLITON SINGLE FOUNDATION EXPORT PROVIDER PACKET GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- n_obligations: `{len(obligations)}`",
                "- Ontology: `NadsolitonSingleFoundation` explicit in both branches.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "n_obligations": len(obligations)}))


if __name__ == "__main__":
    main()
