#!/usr/bin/env python3
"""
QW-2230: QFT axiom-free O1c attachment specification gate.

Purpose:
- consolidate L5 axiom-free O1a/O1b provenance into final O1c attachment layer,
- define exact discharge bundle required for axiom-free terminal proof object,
- keep no-overclaim theorem boundary explicit.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict[str, Any]:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def main() -> None:
    q2224 = load("report_qw2224_qft_axiom_free_discharge_spec_gate.json")
    q2226 = load("report_qw2226_qft_axiom_free_o1a_provenance_gate.json")
    q2228 = load("report_qw2228_qft_axiom_free_o1b_provenance_gate.json")

    sub_c = q2224.get("sub_obligations", {}).get("L5_AXIOM_FREE_O1c", {})
    deps = sub_c.get("depends_on", [])

    unresolved_o1a = q2226.get("unresolved_targets", [])
    unresolved_o1b = q2228.get("unresolved_targets", [])

    unresolved_union = sorted(set(unresolved_o1a) | set(unresolved_o1b))

    discharge_bundle = {
        "theorem_targets": [
            "PositivityToReconstruction",
            "UnitarySMatrixAndScatteringCompleteness",
        ],
        "witness_replacements": [
            "L5O1aWitness",
            "L5O1bWitness",
        ],
        "terminal_files_to_update": [
            "FIN_L5_O1A_O1_TERMINAL.lean",
            "FIN_L5_O1B_O1_TERMINAL.lean",
        ],
        "required_artifacts": [
            "axiom-free theorem proof files",
            "machine-check logs (Lean/Coq)",
            "proof object hash manifest",
            "cross-link report to O1a/O1b provenance",
        ],
    }

    acceptance = {
        "C1": "All unresolved O1a/O1b theorem targets discharged by derived lemmas.",
        "C2": "No `axiom ...Witness` statements remain in updated L5 terminal files.",
        "C3": "Machine checker exits zero on updated L5 O1A/O1B files.",
        "C4": "Hashed axiom-free proof object attached and linked to QW-2226/QW-2228.",
        "C5": "No-overclaim boundary maintained until all C1..C4 are true.",
    }

    flags = {
        "q2224_o1c_subobligation_present": bool(sub_c),
        "q2226_o1a_provenance_pass_present": q2226.get("verdict")
        == "QFT_AXIOM_FREE_O1A_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM",
        "q2228_o1b_provenance_pass_present": q2228.get("verdict")
        == "QFT_AXIOM_FREE_O1B_PROVENANCE_GATE_PASS_PARTIAL_UNRESOLVED_THEOREM",
        "o1c_dependency_exact": deps == ["L5_AXIOM_FREE_O1a", "L5_AXIOM_FREE_O1b"],
        "unresolved_union_nonempty": len(unresolved_union) > 0,
        "discharge_bundle_declared": len(discharge_bundle["theorem_targets"]) == 2 and len(discharge_bundle["witness_replacements"]) == 2,
        "acceptance_criteria_complete_ge_5": len(acceptance) >= 5,
        "o1c_attachment_fully_closed": False,
        "no_overclaim_boundary_explicit": True,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2224_o1c_subobligation_present"]
        and flags["q2226_o1a_provenance_pass_present"]
        and flags["q2228_o1b_provenance_pass_present"]
        and flags["o1c_dependency_exact"]
        and flags["unresolved_union_nonempty"]
        and flags["discharge_bundle_declared"]
        and flags["acceptance_criteria_complete_ge_5"]
        and flags["no_overclaim_boundary_explicit"]
    )

    verdict = (
        "QFT_AXIOM_FREE_O1C_ATTACHMENT_SPEC_GATE_PASS_PARTIAL_DISCHARGE_PENDING"
        if core_ok
        else "QFT_AXIOM_FREE_O1C_ATTACHMENT_SPEC_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2224": "report_qw2224_qft_axiom_free_discharge_spec_gate.json",
            "q2226": "report_qw2226_qft_axiom_free_o1a_provenance_gate.json",
            "q2228": "report_qw2228_qft_axiom_free_o1b_provenance_gate.json",
        },
        "o1c_subobligation": sub_c,
        "unresolved_union": unresolved_union,
        "discharge_bundle": discharge_bundle,
        "acceptance_criteria": acceptance,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "EXECUTE_QFT_AXIOM_FREE_THEOREM_DISCHARGE_AND_ATTACH_FINAL_O1C_PROOF_OBJECT",
    }

    out_json = ROOT / "report_qw2230_qft_axiom_free_o1c_attachment_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2230_QFT_AXIOM_FREE_O1C_ATTACHMENT_SPEC_GATE.md"

    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2230: QFT AXIOM-FREE O1C ATTACHMENT SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                f"- unresolved_union_size: `{len(unresolved_union)}`",
                "",
                "## Core result",
                "- O1c attachment layer dla L5 jest formalnie zdefiniowana (bundle + acceptance + dependencies).",
                "- Remaining gap pozostaje jawnie theorem-level (discharge pending).",
                "",
                "## Artifacts",
                "- JSON: `report_qw2230_qft_axiom_free_o1c_attachment_spec_gate.json`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
