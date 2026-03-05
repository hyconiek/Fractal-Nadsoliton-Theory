#!/usr/bin/env python3
"""
QW-2195: Generation-mapping axiom-augmented gate (L20).

Purpose:
- formalize current L20 boundary after QW-2191..QW-2193:
  * axiom-free physical uniqueness remains open,
  * deterministic axiom-augmented mapping rule is now explicit and auditable.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2195_generation_mapping_axiom_augmented_gate.json"
OUT_MD = ROOT / "RAPORT_QW2195_GENERATION_MAPPING_AXIOM_AUGMENTED_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def overlap_score(labels: List[int], template: Dict[str, List[int]], perm: Tuple[int, int, int]) -> int:
    score = 0
    for tgt in [0, 1, 2]:
        cluster = perm[tgt]
        members = set(template[str(tgt)])
        score += sum(1 for idx, lab in enumerate(labels) if lab == cluster and idx in members)
    return score


def main() -> None:
    r2125 = load_json("report_qw2125_ktotal_generation_alignment_audit.json")
    r2191 = load_json("report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json")
    r2192 = load_json("report_qw2192_mode_index_selection_axiom_gate.json")
    r2193 = load_json("report_qw2193_selection_axiom_family_robustness_gate.json")

    labels = [int(x) for x in r2125["base_tripartition_labels"]]
    tpl_mod3 = r2125["templates"]["mod3_4x3_template"]
    tpl_contig = r2125["templates"]["contiguous_4x3_template"]

    perms = list(itertools.permutations([0, 1, 2]))
    scored = []
    for p in perms:
        s_mod3 = overlap_score(labels, tpl_mod3, p)
        s_contig = overlap_score(labels, tpl_contig, p)
        scored.append({"perm_targets_to_cluster": list(p), "mod3_score_12": s_mod3, "contig_score_12": s_contig})

    best_mod3 = max(scored, key=lambda r: r["mod3_score_12"])["mod3_score_12"]
    best_rows = [r for r in scored if r["mod3_score_12"] == best_mod3]

    # Axiom-augmented deterministic rule:
    # among best-score permutations, choose lexicographically smallest.
    # This is explicit, reproducible, and no-scan over continuous parameters.
    selected = sorted(best_rows, key=lambda r: tuple(r["perm_targets_to_cluster"]))[0]

    flags = {
        "q2125_structural_tripartition_present": bool(r2125["flags"].get("q2118_tripartition_is_balanced_4_4_4", False)),
        "q2125_mod3_template_beats_contiguous": bool(r2125["flags"].get("base_mod3_beats_contiguous_template", False)),
        "q2125_physical_uniqueness_open": bool(not r2125["flags"].get("generation_mapping_is_unique_and_physical", True)),
        "q2191_axiom_free_obstruction_present": bool(str(r2191.get("verdict", "")).endswith("PASS_STRICT")),
        "q2192_axiom_augmented_mode_uniqueness_closed": bool(str(r2192.get("verdict", "")).endswith("UNIQUENESS_CLOSED")),
        "q2193_axiom_family_robustness_present": bool(str(r2193.get("verdict", "")).endswith("AUGMENTED_ROBUST")),
        "axiom_augmented_generation_rule_declared": True,
        "axiom_augmented_generation_mapping_computable": bool(len(best_rows) >= 1),
        "axiom_augmented_generation_mapping_deterministic": True,
        "axiom_augmented_generation_mapping_selected": bool(best_mod3 >= 8),
        "axiom_free_generation_mapping_closed": False,
        "deterministic_no_scan_no_retune": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2125_structural_tripartition_present"]
        and flags["q2125_mod3_template_beats_contiguous"]
        and flags["q2125_physical_uniqueness_open"]
        and flags["q2191_axiom_free_obstruction_present"]
        and flags["q2192_axiom_augmented_mode_uniqueness_closed"]
        and flags["q2193_axiom_family_robustness_present"]
        and flags["axiom_augmented_generation_rule_declared"]
        and flags["axiom_augmented_generation_mapping_computable"]
        and flags["axiom_augmented_generation_mapping_deterministic"]
        and flags["axiom_augmented_generation_mapping_selected"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "GENERATION_MAPPING_AXIOM_AUGMENTED_GATE_PASS_PARTIAL_AXIOM_FREE_OPEN"
        if core_ok
        else "GENERATION_MAPPING_AXIOM_AUGMENTED_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2125": "report_qw2125_ktotal_generation_alignment_audit.json",
            "q2191": "report_qw2191_mode_index_uniqueness_obstruction_theorem_gate.json",
            "q2192": "report_qw2192_mode_index_selection_axiom_gate.json",
            "q2193": "report_qw2193_selection_axiom_family_robustness_gate.json",
        },
        "axiom_augmented_rule": {
            "name": "max_mod3_overlap_then_lexicographic_tie_break",
            "description": "Choose permutation with maximal mod3 overlap score; if multiple, choose lexicographically smallest.",
            "best_mod3_score_12": best_mod3,
            "n_best_permutations": len(best_rows),
            "selected_permutation_targets_to_cluster": selected["perm_targets_to_cluster"],
        },
        "audit": {
            "base_labels": labels,
            "permutation_scores": scored,
            "q2125_perturbation_robustness": r2125.get("perturbation_robustness", {}),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROVIDE_PHYSICAL_JUSTIFICATION_FOR_AXIOM_AUGMENTED_GENERATION_RULE_OR_DERIVE_AXIOM_FREE_MAPPING"
            if verdict.endswith("AXIOM_FREE_OPEN")
            else "REPAIR_GENERATION_MAPPING_AXIOM_AUGMENTED_CHAIN_AND_RERUN_QW2195"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2195: GENERATION MAPPING AXIOM AUGMENTED GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- L20 remains open in axiom-free physical uniqueness scope (explicitly).",
        "- A deterministic axiom-augmented generation mapping rule is now explicit and auditable.",
        "- Rule is chained to QW-2191..2193 boundaries and robustness checks.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
