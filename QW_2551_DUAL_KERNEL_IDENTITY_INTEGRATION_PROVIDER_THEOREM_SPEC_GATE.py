#!/usr/bin/env python3
"""QW-2551: dual kernel-identity-integration provider theorem spec gate."""

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
    q2550 = load("report_qw2550_dual_kernel_identity_integration_provider_minimal_blocker_cut_gate.json")
    q2549 = load("report_qw2549_strict_anti_false_pass_identity_integration_frontier_gate.json")

    acceptance_criteria = {
        "C1": "Identity-integration domain invariance theorem in declared operator class.",
        "C2": "Self-adjointness/positivity preservation under integration constraints.",
        "C3": "Integration coercive lower-bound theorem with explicit constants.",
        "C4": "Spectral/graph persistence under bounded integration-preserving perturbations.",
        "C5": "Bridge theorem: integration lemmas imply well-posedness/positivity targets.",
        "C6": "Machine-check package with reproducible proof hashes and dependency DAG.",
    }

    theorem_targets = {
        "L12_target": "RG_KernelIdentityIntegrationToWellPosedness_Theorem",
        "L5_target": "QFT_KernelIdentityIntegrationToPositivity_Theorem",
    }

    lemma_dag = {
        "RG_KIINT_L1_IntegrationDomainInvariance": [],
        "RG_KIINT_L2_SelfAdjointPositivityPreservation": ["RG_KIINT_L1_IntegrationDomainInvariance"],
        "RG_KIINT_L3_IntegrationCoerciveLowerBound": ["RG_KIINT_L2_SelfAdjointPositivityPreservation"],
        "RG_KIINT_L4_BoundedIntegrationPreservation": ["RG_KIINT_L3_IntegrationCoerciveLowerBound"],
        "RG_KernelIdentityIntegrationToWellPosedness_Theorem": ["RG_KIINT_L4_BoundedIntegrationPreservation"],
        "QFT_KIINT_L1_IntegrationDomainInvariance": [],
        "QFT_KIINT_L2_SelfAdjointPositivityPreservation": ["QFT_KIINT_L1_IntegrationDomainInvariance"],
        "QFT_KIINT_L3_IntegrationCoerciveLowerBound": ["QFT_KIINT_L2_SelfAdjointPositivityPreservation"],
        "QFT_KIINT_L4_BoundedIntegrationPreservation": ["QFT_KIINT_L3_IntegrationCoerciveLowerBound"],
        "QFT_KernelIdentityIntegrationToPositivity_Theorem": ["QFT_KIINT_L4_BoundedIntegrationPreservation"],
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
            "statement": "Integration-constrained domain is dense and graph-closed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_2": {
            "class": "technical",
            "statement": "Self-adjoint extension uniqueness in declared integration class.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_3": {
            "class": "technical",
            "statement": "Explicit coercive lower-bound constants in integration norm.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_4": {
            "class": "technical",
            "statement": "Perturbation bound below integration-preservation threshold.",
            "applies_to": ["L12_target", "L5_target"],
        },
    }

    rg_nodes = [node for node in lemma_dag if node.startswith("RG_")]
    qft_nodes = [node for node in lemma_dag if node.startswith("QFT_")]
    acyclic = dag_is_acyclic(lemma_dag)

    flags = {
        "q2550_min_cut_ready": q2550.get("verdict")
        == "DUAL_KERNEL_IDENTITY_INTEGRATION_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2549_integration_ready": q2549.get("verdict")
        == "STRICT_ANTI_FALSE_PASS_IDENTITY_INTEGRATION_FRONTIER_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
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
        flags["q2550_min_cut_ready"]
        and flags["q2549_integration_ready"]
        and flags["theorem_targets_count_is_two"]
        and flags["acceptance_criteria_complete_ge_6"]
        and flags["lemma_dag_acyclic"]
        and flags["rg_lemma_nodes_ge_5"]
        and flags["qft_lemma_nodes_ge_5"]
        and flags["assumption_map_complete_ge_6"]
        and flags["assumption_classification_has_physical_and_technical"]
    )

    verdict = (
        "DUAL_KERNEL_IDENTITY_INTEGRATION_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "DUAL_KERNEL_IDENTITY_INTEGRATION_PROVIDER_THEOREM_SPEC_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2550": "report_qw2550_dual_kernel_identity_integration_provider_minimal_blocker_cut_gate.json",
            "q2549": "report_qw2549_strict_anti_false_pass_identity_integration_frontier_gate.json",
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
    spec_path = ROOT / "spec_qw2551_dual_kernel_identity_integration_provider_theorem_spec.json"
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
    proof_path = ROOT / "proof_object_qw2551_dual_kernel_identity_integration_provider_theorem_spec.json"
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
        "required_next_step": "RUN_DUAL_KERNEL_IDENTITY_INTEGRATION_PROVIDER_COUNTEREXAMPLE_SEARCH",
    }

    out_json = ROOT / "report_qw2551_dual_kernel_identity_integration_provider_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2551_DUAL_KERNEL_IDENTITY_INTEGRATION_PROVIDER_THEOREM_SPEC_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2551: DUAL KERNEL IDENTITY INTEGRATION PROVIDER THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zdefiniowano 2 target twierdzenia aktualnego frontu identity-integration provider.",
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
