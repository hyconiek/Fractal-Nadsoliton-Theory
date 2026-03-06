#!/usr/bin/env python3
"""QW-2461: dual kernel-spectral-invariance provider theorem spec gate.

Builds explicit theorem specification for current dual frontier:
- RG_KernelSpectralInvarianceToWellPosedness_Theorem
- QFT_KernelSpectralInvarianceToPositivity_Theorem

Includes minimal lemma decomposition and explicit assumption map
with physical vs technical classification.
No theorem-level/full-closure claim is allowed.
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

    for n in graph:
        if n not in perm_mark and not visit(n):
            return False
    return True


def main() -> None:
    q2459 = load("report_qw2459_dual_kernel_spectral_invariance_provider_minimal_blocker_cut_gate.json")
    q2460 = load("report_qw2460_strict_anti_false_pass_spectral_chain_continuation_gate.json")

    acceptance_criteria = {
        "C1": "Domain invariance theorem for spectral flow operator in declared scope.",
        "C2": "Self-adjointness preservation on invariant domain under declared transforms.",
        "C3": "Positivity/coercivity invariance bound with explicit constants.",
        "C4": "Spectral stability under bounded perturbation in invariance regime.",
        "C5": "Bridge theorem from invariance-lemmas to well-posedness/positivity targets.",
        "C6": "Machine-check package with reproducible proof-object hashes and dependency DAG.",
    }

    theorem_targets = {
        "L12_target": "RG_KernelSpectralInvarianceToWellPosedness_Theorem",
        "L5_target": "QFT_KernelSpectralInvarianceToPositivity_Theorem",
    }

    lemma_dag = {
        "RG_KSI_L1_DomainInvariance": [],
        "RG_KSI_L2_SelfAdjointPreservation": ["RG_KSI_L1_DomainInvariance"],
        "RG_KSI_L3_CoerciveInvarianceBound": ["RG_KSI_L2_SelfAdjointPreservation"],
        "RG_KSI_L4_PerturbativeSpectralStability": ["RG_KSI_L3_CoerciveInvarianceBound"],
        "RG_KernelSpectralInvarianceToWellPosedness_Theorem": [
            "RG_KSI_L4_PerturbativeSpectralStability"
        ],
        "QFT_KSI_L1_DomainInvariance": [],
        "QFT_KSI_L2_SelfAdjointPreservation": ["QFT_KSI_L1_DomainInvariance"],
        "QFT_KSI_L3_PositivityInvarianceBound": ["QFT_KSI_L2_SelfAdjointPreservation"],
        "QFT_KSI_L4_PerturbativePositivityStability": ["QFT_KSI_L3_PositivityInvarianceBound"],
        "QFT_KernelSpectralInvarianceToPositivity_Theorem": [
            "QFT_KSI_L4_PerturbativePositivityStability"
        ],
    }

    assumptions = {
        "A_phys_1": {
            "class": "physical",
            "statement": "Complete FIN action and declared branch are fixed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_phys_2": {
            "class": "physical",
            "statement": "Constructive RG/QFT map exists in declared strict scope.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_1": {
            "class": "technical",
            "statement": "Operator domain is dense and graph-closed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_2": {
            "class": "technical",
            "statement": "Self-adjoint extension exists and is unique in declared class.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_3": {
            "class": "technical",
            "statement": "Lower-bound positivity/coercivity constant is explicit and strict.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_4": {
            "class": "technical",
            "statement": "Perturbation norm bound remains below invariance threshold.",
            "applies_to": ["L12_target", "L5_target"],
        },
    }

    rg_nodes = [n for n in lemma_dag if n.startswith("RG_")]
    qft_nodes = [n for n in lemma_dag if n.startswith("QFT_")]
    acyclic = dag_is_acyclic(lemma_dag)

    flags = {
        "q2459_min_cut_ready": q2459.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
        "q2460_integrity_ready": q2460.get("verdict")
        == "STRICT_ANTI_FALSE_PASS_SPECTRAL_CHAIN_CONTINUATION_GATE_PASS_WITH_BLOCKERS_EXPLICIT",
        "theorem_targets_count_is_two": len(theorem_targets) == 2,
        "acceptance_criteria_complete_ge_6": len(acceptance_criteria) >= 6,
        "lemma_dag_acyclic": acyclic,
        "rg_lemma_nodes_ge_5": len(rg_nodes) >= 5,
        "qft_lemma_nodes_ge_5": len(qft_nodes) >= 5,
        "assumption_map_complete_ge_6": len(assumptions) >= 6,
        "assumption_classification_has_physical_and_technical": {
            assumptions[k]["class"] for k in assumptions
        }
        == {"physical", "technical"},
        "all_strict_obligations_fully_closed": False,
    }

    pass_count = sum(1 for v in flags.values() if v)
    total_flags = len(flags)

    core_ok = (
        flags["q2459_min_cut_ready"]
        and flags["q2460_integrity_ready"]
        and flags["theorem_targets_count_is_two"]
        and flags["acceptance_criteria_complete_ge_6"]
        and flags["lemma_dag_acyclic"]
        and flags["rg_lemma_nodes_ge_5"]
        and flags["qft_lemma_nodes_ge_5"]
        and flags["assumption_map_complete_ge_6"]
        and flags["assumption_classification_has_physical_and_technical"]
    )

    verdict = (
        "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREM_SPEC_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2459": "report_qw2459_dual_kernel_spectral_invariance_provider_minimal_blocker_cut_gate.json",
            "q2460": "report_qw2460_strict_anti_false_pass_spectral_chain_continuation_gate.json",
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
    spec_path = ROOT / "spec_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec.json"
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
    proof_path = ROOT / "proof_object_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec.json"
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
        "required_next_step": "RUN_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_COUNTEREXAMPLE_SEARCH",
    }

    out_json = ROOT / "report_qw2461_dual_kernel_spectral_invariance_provider_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2461_DUAL_KERNEL_SPECTRAL_INVARIANCE_PROVIDER_THEOREM_SPEC_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2461: DUAL KERNEL SPECTRAL INVARIANCE PROVIDER THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zdefiniowano 2 target twierdzenia aktualnego frontu invariance-provider.",
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
