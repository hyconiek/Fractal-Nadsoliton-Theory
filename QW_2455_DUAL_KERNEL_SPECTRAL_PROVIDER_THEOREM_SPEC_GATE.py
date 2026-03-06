#!/usr/bin/env python3
"""QW-2455: dual kernel-spectral provider theorem spec gate.

Builds explicit theorem specification for current dual frontier:
- RG_KernelSpectralClosureToWellPosedness_Theorem
- QFT_KernelSpectralClosureToPositivity_Theorem

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
    q2453 = load("report_qw2453_non_axiomatic_dual_kernel_operator_closure_provider_derivation_attempt_gate.json")
    q2454 = load("report_qw2454_dual_kernel_spectral_provider_minimal_blocker_cut_gate.json")

    acceptance_criteria = {
        "C1": "Self-adjoint closure of effective spectral operator in declared domain.",
        "C2": "Strict coercive lower-bound theorem with explicit alpha>0 constants.",
        "C3": "Bounded perturbation theorem preserving positivity/well-posedness under explicit norm control.",
        "C4": "Domain-closure and graph-norm continuity theorem for operator extension.",
        "C5": "Composition theorem: spectral closure lemmas imply RG well-posedness and QFT positivity targets.",
        "C6": "Machine-check package with reproducible proof-object hashes and dependency DAG.",
    }

    theorem_targets = {
        "L12_target": "RG_KernelSpectralClosureToWellPosedness_Theorem",
        "L5_target": "QFT_KernelSpectralClosureToPositivity_Theorem",
    }

    lemma_dag = {
        "RG_KS_L1_SelfAdjointClosure": [],
        "RG_KS_L2_CoerciveLowerBound": ["RG_KS_L1_SelfAdjointClosure"],
        "RG_KS_L3_BoundedPerturbationControl": ["RG_KS_L2_CoerciveLowerBound"],
        "RG_KS_L4_DomainClosureContinuity": ["RG_KS_L1_SelfAdjointClosure"],
        "RG_KernelSpectralClosureToWellPosedness_Theorem": [
            "RG_KS_L3_BoundedPerturbationControl",
            "RG_KS_L4_DomainClosureContinuity",
        ],
        "QFT_KS_L1_SelfAdjointClosure": [],
        "QFT_KS_L2_PositiveLowerBound": ["QFT_KS_L1_SelfAdjointClosure"],
        "QFT_KS_L3_BoundedPerturbationPositivity": ["QFT_KS_L2_PositiveLowerBound"],
        "QFT_KS_L4_DomainClosureContinuity": ["QFT_KS_L1_SelfAdjointClosure"],
        "QFT_KernelSpectralClosureToPositivity_Theorem": [
            "QFT_KS_L3_BoundedPerturbationPositivity",
            "QFT_KS_L4_DomainClosureContinuity",
        ],
    }

    assumptions = {
        "A_phys_1": {
            "class": "physical",
            "statement": "Complete FIN action and declared branch selection are fixed.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_phys_2": {
            "class": "physical",
            "statement": "Constructive map/scheme exists in declared strict scope.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_1": {
            "class": "technical",
            "statement": "Operators are densely-defined and symmetric/self-adjoint on stated domain.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_2": {
            "class": "technical",
            "statement": "Strict lower-bound alpha>0 is explicit for unperturbed operator.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_3": {
            "class": "technical",
            "statement": "Perturbation norm bound is explicit and below stability threshold.",
            "applies_to": ["L12_target", "L5_target"],
        },
        "A_tech_4": {
            "class": "technical",
            "statement": "Domain closure continuity under graph norm is proven.",
            "applies_to": ["L12_target", "L5_target"],
        },
    }

    rg_nodes = [n for n in lemma_dag if n.startswith("RG_")]
    qft_nodes = [n for n in lemma_dag if n.startswith("QFT_")]
    acyclic = dag_is_acyclic(lemma_dag)

    flags = {
        "q2453_attempt_blocked_at_kernel_spectral_layer": q2453.get("verdict")
        == "NON_AXIOMATIC_DUAL_KERNEL_OPERATOR_CLOSURE_PROVIDER_DERIVATION_ATTEMPT_GATE_PASS_PARTIAL_BLOCKED_BY_KERNEL_SPECTRAL_CLOSURE_PROVIDER_THEOREMS",
        "q2454_min_cut_ready": q2454.get("verdict")
        == "DUAL_KERNEL_SPECTRAL_PROVIDER_MINIMAL_BLOCKER_CUT_GATE_PASS_PARTIAL_TWO_SYMBOLS_ISOLATED",
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
        flags["q2453_attempt_blocked_at_kernel_spectral_layer"]
        and flags["q2454_min_cut_ready"]
        and flags["theorem_targets_count_is_two"]
        and flags["acceptance_criteria_complete_ge_6"]
        and flags["lemma_dag_acyclic"]
        and flags["rg_lemma_nodes_ge_5"]
        and flags["qft_lemma_nodes_ge_5"]
        and flags["assumption_map_complete_ge_6"]
        and flags["assumption_classification_has_physical_and_technical"]
    )

    verdict = (
        "DUAL_KERNEL_SPECTRAL_PROVIDER_THEOREM_SPEC_GATE_PASS_PARTIAL_TERMINAL_LAYER_READY"
        if core_ok
        else "DUAL_KERNEL_SPECTRAL_PROVIDER_THEOREM_SPEC_GATE_FAIL"
    )

    spec = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "report_qw2454_dual_kernel_spectral_provider_minimal_blocker_cut_gate.json",
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
    spec_path = ROOT / "spec_qw2455_dual_kernel_spectral_provider_theorem_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    proof_obj = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2453": "report_qw2453_non_axiomatic_dual_kernel_operator_closure_provider_derivation_attempt_gate.json",
            "q2454": "report_qw2454_dual_kernel_spectral_provider_minimal_blocker_cut_gate.json",
        },
        "counts": {
            "n_acceptance_criteria": len(acceptance_criteria),
            "n_lemma_nodes_total": len(lemma_dag),
            "n_assumptions": len(assumptions),
        },
        "lemma_dag_acyclic": acyclic,
    }
    proof_path = ROOT / "proof_object_qw2455_dual_kernel_spectral_provider_theorem_spec.json"
    proof_path.write_text(json.dumps(proof_obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": proof_obj["sources"],
        "spec_file": spec_path.name,
        "spec_sha256": sha256_file(spec_path),
        "proof_object_file": proof_path.name,
        "proof_object_sha256": sha256_file(proof_path),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": "RUN_DUAL_KERNEL_SPECTRAL_PROVIDER_COUNTEREXAMPLE_SEARCH",
    }

    out_json = ROOT / "report_qw2455_dual_kernel_spectral_provider_theorem_spec_gate.json"
    out_md = ROOT / "RAPORT_QW2455_DUAL_KERNEL_SPECTRAL_PROVIDER_THEOREM_SPEC_GATE.md"
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    out_md.write_text(
        "\n".join(
            [
                "# RAPORT QW-2455: DUAL KERNEL SPECTRAL PROVIDER THEOREM SPEC GATE",
                "",
                f"- Date UTC: {out['generated_utc']}",
                f"- Verdict: **{verdict}**",
                f"- pass_count: `{pass_count}/{total_flags}`",
                "",
                "## Core result",
                "- Zdefiniowano 2 target twierdzen na aktualnym froncie blockerow.",
                "- Zbudowano minimalny DAG lematow (RG + QFT) oraz jawna mape zalozen (physical/technical).",
                "- Brak podstaw do theorem-level/full-closure PASS.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}))


if __name__ == "__main__":
    main()
