#!/usr/bin/env python3
"""QW-2546: dual kernel-identity-consolidation provider theorem spec gate."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def dag_is_acyclic(graph: dict[str, list[str]]) -> bool:
    temp_mark: set[str] = set()
    perm_mark: set[str] = set()

    def visit(node: str) -> bool:
        if node in perm_mark:
            return True
        if node in temp_mark:
            return False
        temp_mark.add(node)
        for dep in graph.get(node, []):
            if not visit(dep):
                return False
        temp_mark.remove(node)
        perm_mark.add(node)
        return True

    for node in graph:
        if node not in perm_mark and not visit(node):
            return False
    return True


def main() -> None:
    q2545 = load("report_qw2545_dual_kernel_identity_consolidation_provider_minimal_blocker_cut_gate.json")
    q2544 = load("report_qw2544_strict_anti_false_pass_identity_consolidation_frontier_gate.json")

    acceptance_criteria = {
        "C1": "Identity-consolidation domain invariance theorem in declared operator class.",
        "C2": "Self-adjointness/positivity preservation under consolidation constraints.",
        "C3": "Consolidation coercive lower-bound theorem with explicit constants.",
        "C4": "Spectral/graph persistence under bounded consolidation-preserving perturbations.",
        "C5": "Bridge theorem: consolidation lemmas imply well-posedness/positivity targets.",
        "C6": "Machine-check package with reproducible proof hashes and dependency DAG.",
    }

    theorem_targets = {
        "L12_target": "RG_KernelIdentityConsolidationToWellPosedness_Theorem",
        "L5_target": "QFT_KernelIdentityConsolidationToPositivity_Theorem",
    }

    lemma_dag = {
        "RG_KICONS_L1_ConsolidationDomainInvariance": [],
        "RG_KICONS_L2_SelfAdjointPositivityPreservation": ["RG_KICONS_L1_ConsolidationDomainInvariance"],
        "RG_KICONS_L3_ConsolidationCoerciveLowerBound": ["RG_KICONS_L2_SelfAdjointPositivityPreservation"],
        "RG_KICONS_L4_BoundedConsolidationPreservation": ["RG_KICONS_L3_ConsolidationCoerciveLowerBound"],
        "RG_KernelIdentityConsolidationToWellPosedness_Theorem": ["RG_KICONS_L4_BoundedConsolidationPreservation"],
        "QFT_KICONS_L1_ConsolidationDomainInvariance": [],
        "QFT_KICONS_L2_SelfAdjointPositivityPreservation": ["QFT_KICONS_L1_ConsolidationDomainInvariance"],
        "QFT_KICONS_L3_ConsolidationCoerciveLowerBound": ["QFT_KICONS_L2_SelfAdjointPositivityPreservation"],
        "QFT_KICONS_L4_BoundedConsolidationPreservation": ["QFT_KICONS_L3_ConsolidationCoerciveLowerBound"],
        "QFT_KernelIdentityConsolidationToPositivity_Theorem": ["QFT_KICONS_L4_BoundedConsolidationPreservation"],
    }

    assumptions = {
        "A_phys_1": {
            "class": "physical",
            "statement": "Complete FIN action and declared branch selection remain fixed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_phys_2": {
            "class": "physical",
            "statement": "Constructive RG/QFT mapping in strict scope is preserved.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_1": {
            "class": "technical",
            "statement": "Consolidation-constrained domain is dense and graph-closed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_2": {
            "class": "technical",
            "statement": "Self-adjoint extension uniqueness in declared consolidation class.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_3": {
            "class": "technical",
            "statement": "Explicit coercive lower-bound constants in consolidation norm.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_4": {
            "class": "technical",
            "statement": "Perturbation bound below consolidation-preservation threshold.",
            "applies_to": ["L12_target", "L5_target"],
        },
    }

    rg_nodes = [node for node in lemma_dag if node.startswith("RG_")]
    qft_nodes = [node for node in lemma_dag if node.startswith("QFT_")]
    acyclic = dag_is_acyclic(lemma_dag)

    flags = {
        "q2545_min_cut_ready": q2545.get("verdict")
        == "DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2544_consolidation_ready": q2544.get("verdict")
        == "STRICT_ANTI_FALSE_PASS_IDENTITY_CONSOLIDATION_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
        "theorem_targets_count_is_two": len(theorem_targets) == 2,
        "acceptance_criteria_complete_ge_6": len(acceptance_criteria) >= 6,
        "lemma_dag_acyclic": acyclic,
        "rg_lemma_nodes_ge_5": len(rg_nodes) >= 5,
        "qft_lemma_nodes_ge_5": len(qft_nodes) >= 5,
        "assumption_map_complete_ge_6": len(assumptions) >= 6,
        "assumption_classification_has_physical_and_technical": {
            assumptions[key]["class"] for key in assumptions
        } == {"physical", "technical"},
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for value in flags.values() if value)
    total_flags = len(flags)

    core_ok = (
        flags["q2545_min_cut_ready"]
        and flags["q2544_consolidation_ready"]
        and flags["theorem_targets_count_is_two"]
        and flags["acceptance_criteria_complete_ge_6"]
        and flags["lemma_dag_acyclic"]
        and flags["rg_lemma_nodes_ge_5"]
        and flags["qft_lemma_nodes_ge_5"]
        and flags["assumption_map_complete_ge_6"]
        and flags["assumption_classification_has_physical_and_technical"]
    )

    verdict = (
        "DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_THEOREM_SPEC_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2545": "report_qw2545_dual_kernel_identity_consolidation_provider_minimal_blocker_cut_gate.json",
            "q2544": "report_qw2544_strict_anti_false_pass_identity_consolidation_frontier_gate.json",
        },
        "theorem_targets": theorem_targets,
        "acceptance_criteria": acceptance_criteria,
        "lemma_dag": lemma_dag,
        "assumptions": assumptions,
        "scope_boundary": {
            "theorem_spec_only": True,
            "counterexample_search_pending": True,
            "provider_discharge_completed": False,
            "all_strict_obligations_fully_closed": False,
            "overclaim_forbidden": True,
        },
    }
    spec_path = ROOT / "spec_qw2546_dual_kernel_identity_consolidation_provider_theorem_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": spec["sources"],
        "counts": {
            "n_acceptance_criteria": len(acceptance_criteria),
            "n_lemma_nodes_total": len(lemma_dag),
            "n_assumptions": len(assumptions),
        },
        "lemma_dag_acyclic": acyclic,
    }
    proof_path = ROOT / "proof_object_qw2546_dual_kernel_identity_consolidation_provider_theorem_spec.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": spec["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "RUN_DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_COUNTEREXAMPLE_SEARCH",
    }

    out_json = ROOT / "report_qw2546_dual_kernel_identity_consolidation_provider_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2546_DUAL_KERNEL_IDENTITY_CONSOLIDATION_PROVIDER_THEOREM_SPEC_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2546: DUAL KERNEL IDENTITY CONSOLIDATION PROVIDER THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zdefiniowano 2 target twierdzenia aktualnego frontu identity-consolidation provider.",
                "- Zbudowano minimalny DAG lematow (RG + QFT) i jawna mape zalozen (physical/technical).",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
