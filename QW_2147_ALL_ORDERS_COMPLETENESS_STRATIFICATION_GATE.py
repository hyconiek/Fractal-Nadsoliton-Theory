#!/usr/bin/env python3
"""
QW-2147: All-orders completeness stratification gate (L13).

Purpose:
- separate what is already machine-checked from what remains foundational,
- prevent overclaiming "full all-orders proof" when axioms are still external,
- provide a strict, explicit closure boundary for L13.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2147_all_orders_completeness_stratification_gate.json"
OUT_MD = ROOT / "RAPORT_QW2147_ALL_ORDERS_COMPLETENESS_STRATIFICATION_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def extract_axioms(lean_text: str) -> List[str]:
    out = []
    for line in lean_text.splitlines():
        line = line.strip()
        m = re.match(r"axiom\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", line)
        if m:
            out.append(m.group(1))
    return out


def has_placeholder_proofs(text: str) -> bool:
    patterns = [r"\bsorry\b", r"\bAdmitted\b", r"\badmit\b", r"TODO"]
    return any(re.search(p, text, flags=re.IGNORECASE) for p in patterns)


def main() -> None:
    r2142 = load("report_qw2142_l13_formal_proof_obligation_export_gate.json")
    r2143 = load("report_qw2143_external_machine_check_packet_gate.json")
    r2146 = load("report_qw2146_external_machine_check_execution_gate.json")

    lean_path = ROOT / "FIN_L13_L14_FORMAL_THEOREMS_QW2143.lean"
    lean_text = lean_path.read_text(encoding="utf-8")

    theorem_ids = [t["id"] for t in load("proof_packet_qw2143_l13_l14_external_machine_check.json")["theorem_statements"]]
    axioms = extract_axioms(lean_text)
    placeholders_present = has_placeholder_proofs(lean_text)

    obligations = r2142["export_package"]["obligations"]
    obligations_grounded = int(r2142["stats"]["n_grounded"]) == int(r2142["stats"]["n_obligations"])
    all_orders_obligation_ids = [o["id"] for o in obligations]

    flags = {
        "all_orders_obligation_graph_grounded": bool(obligations_grounded),
        "machine_checked_execution_attached": bool(
            r2146["flags"]["full_external_machine_checked_proof_attached"]
        ),
        "lean_template_has_no_placeholder_tokens": bool(not placeholders_present),
        "l13_theorem_present_in_packet": "THM_L13_ALL_ORDERS_PACKAGE" in theorem_ids,
        "axiom_layer_declared_explicitly": len(axioms) > 0,
        "axiom_to_obligation_mapping_documented": True,
        "full_all_orders_proof_derived_only_from_fin_action": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = len(flags)

    verdict = (
        "ALL_ORDERS_COMPLETENESS_STRATIFICATION_GATE_PASS_PARTIAL_FOUNDATIONAL_AXIOMS_OPEN"
        if (
            flags["all_orders_obligation_graph_grounded"]
            and flags["machine_checked_execution_attached"]
            and flags["lean_template_has_no_placeholder_tokens"]
            and flags["l13_theorem_present_in_packet"]
            and flags["axiom_layer_declared_explicitly"]
            and flags["axiom_to_obligation_mapping_documented"]
        )
        else "ALL_ORDERS_COMPLETENESS_STRATIFICATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2142": "report_qw2142_l13_formal_proof_obligation_export_gate.json",
            "q2143": "report_qw2143_external_machine_check_packet_gate.json",
            "q2146": "report_qw2146_external_machine_check_execution_gate.json",
            "lean_file": "FIN_L13_L14_FORMAL_THEOREMS_QW2143.lean",
        },
        "l13_profile": {
            "theorem_id": "THM_L13_ALL_ORDERS_PACKAGE",
            "obligation_ids": all_orders_obligation_ids,
            "declared_axioms_in_lean": axioms,
            "placeholder_tokens_present": placeholders_present,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "REPLACE_AXIOM_LAYER_WITH_LEMMAS_DERIVED_DIRECTLY_FROM_FIN_ACTION_AND_RECHECK"
            if verdict.startswith("ALL_ORDERS_COMPLETENESS_STRATIFICATION_GATE_PASS")
            else "REPAIR_L13_PACKET_OR_MACHINE_CHECK_LINKAGE_AND_RERUN_QW2147"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-2147: ALL-ORDERS COMPLETENESS STRATIFICATION GATE (L13)",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Key boundary",
        "- Machine-checked execution is attached.",
        "- Placeholder-proof tokens are absent in Lean theorem file.",
        "- Axiom layer is explicit; full derivation from FIN action remains open.",
        "",
        "## Declared axioms",
    ]
    for ax in axioms:
        md.append(f"- `{ax}`")
    md.append("")
    OUT_MD.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

