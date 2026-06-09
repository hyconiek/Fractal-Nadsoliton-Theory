#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.json"
MD = GEN / "p2636_s1586_current_toe_blocker_lattice_full_kernel_decision_audit.md"
REPO_ROOT = REPO

SOURCE_FILES = {
    "P2625_NONLINEAR_DAMPING_COMPLETION_SOURCE": GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json",
    "P2629_ZBETA_NORMALIZATION_GAUGE": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2630_BETA_SOURCE_SEPARATION": GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json",
    "P2631_NEURAL_BETA_CRITICALITY": GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json",
    "P2633_DIAGRAM_GROUNDED_RETENTION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
    "P2634_STABILITY_VS_ROLE_COMPLETENESS": GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json",
    "P2635_TOE_NEURAL_UNIVERSE_SIGNATURE": GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json",
    "QW2191_FORMAL_CHAIN": GEN / "p1294_qw2191_r10_formal_proof_chain_and_countermodel_sweep_checkpoint.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "full_strict_kernel_claimed",
    "toe_closure_claimed",
    "legacy_to_strict_bridge_completion_exported",
    "legacy_role_transfer_revalidated",
    "selector_source_exported",
    "q_w_2191_discharged",
    "positive_beta_renormalization_source_exported",
    "phase_frequency_node_gauge_certificate_exported",
    "inverse_hierarchy_role_transfer_exported",
    "role_bearing_ltotal_reenabled",
]


def sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:70]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: full-kernel and closure blockers by research meaning, not packet IDs.
    patterns = {
        "full_kernel_finality_content": (
            "full kernel|complete kernel|final.*kernel|kernel finality|strict kernel.*final|"
            "pełny kernel|domknięcie teorii|ToE closure|single kernel.*ToE"
        ),
        "bridge_completion_content": (
            "legacy.*strict.*bridge|bridge completion|completion bridge|completion certificate|"
            "role-transfer audit|characteristic-by-characteristic|K_legacy_ont|K_strict_gate"
        ),
        "beta_source_content": (
            "positive_beta_renormalization_source|beta normalization|Z_beta|UV normalization|"
            "information conservation beta|edge of chaos|beta_tors"
        ),
        "phase_node_inverse_hierarchy_content": (
            "phase-frequency|node/gauge|integer node|zero lattice|inverse hierarchy|"
            "distant octave|Wilson-loop|cos\\(omega|cos\\(ω"
        ),
        "selector_source_qw2191_content": (
            "QW-2191|selector source|strict-core selector|symmetry-breaking|orientation source|"
            "source-side quotient|global selector"
        ),
        "empirical_neural_universe_content": (
            "neural universe|samoucząca|heavy-tailed attention|positional encoding|CMB/LSS|"
            "GW/PTA|frozen-kernel|holdout|preregistration|external confirmatory"
        ),
        "ltotal_dynamics_content": (
            "role-bearing L_total|L_total|dynamical source|variational functional|"
            "minimum-roughness|Euler-Lagrange|RG stationarity"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for current ToE blockers and full-kernel decision",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def source_fingerprints() -> dict[str, Any]:
    return {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()}


def bool_get(data: dict[str, Any], path: list[str], default: bool = False) -> bool:
    cur: Any = data
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return bool(cur)


def current_blocker_lattice(src: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p2625 = src["P2625_NONLINEAR_DAMPING_COMPLETION_SOURCE"]
    p2629 = src["P2629_ZBETA_NORMALIZATION_GAUGE"]
    p2630 = src["P2630_BETA_SOURCE_SEPARATION"]
    p2631 = src["P2631_NEURAL_BETA_CRITICALITY"]
    p2633 = src["P2633_DIAGRAM_GROUNDED_RETENTION"]
    p2634 = src["P2634_STABILITY_VS_ROLE_COMPLETENESS"]
    p2635 = src["P2635_TOE_NEURAL_UNIVERSE_SIGNATURE"]

    p2633_failed = p2633.get("source_acceptance", {}).get("failed_gates", [])
    p2634_scores = p2634.get("aggregate_scores", {})
    p2635_acceptance = p2635.get("source_acceptance", {})
    p2635_failed = p2635_acceptance.get("failed_gates", [])

    blockers = [
        {
            "blocker": "legacy_to_strict_bridge_completion",
            "why_it_blocks_full_kernel": "A full strict kernel must be licensed as the completed/enriched continuation of the legacy bridge kernel, not only as a stable independent ansatz.",
            "current_evidence": "P2633/P2634 retain strict as an enriched working successor but not a complete role-preserving bridge.",
            "is_closed": False,
            "severity_0_to_1": 1.0,
            "evidence_signals": {
                "p2634_open_role_blocks": p2634_scores.get("role_completeness_open_blocks"),
                "p2633_accepts_all_diagram_characteristics": p2633.get("source_acceptance", {}).get("accepts_all_diagram_characteristics_preserved_by_strict"),
            },
            "next_discharge_test": "Build a characteristic-by-characteristic completion certificate with explicit residual strict-side additions.",
        },
        {
            "blocker": "positive_beta_or_uv_normalization_source",
            "why_it_blocks_full_kernel": "The strict damping coefficient and bridge Z_beta cannot be promoted by convention, internal beta evidence, or neural edge-of-chaos language alone.",
            "current_evidence": "P2629/P2630/P2631 reject target-independent beta/Z_beta source closure.",
            "is_closed": False,
            "severity_0_to_1": 0.95,
            "evidence_signals": {
                "p2629_status": p2629.get("status"),
                "p2630_bridge_source_truth_table_rows": len(p2630.get("bridge_source_truth_table", [])),
                "p2631_information_conservation_beta_identity": p2631.get("source_acceptance", {}).get("accepts_information_conservation_beta_identity"),
            },
            "next_discharge_test": "Derive a dimensionless conservation/normalization identity before comparing to beta=1 or Z_beta=100.",
        },
        {
            "blocker": "phase_frequency_node_gauge_certificate",
            "why_it_blocks_full_kernel": "The phase channel can be stable and useful while still failing the exact legacy node/gauge interpretation.",
            "current_evidence": "P2633 found the declared integer nodes are not exact zeros of cos(pi*d/4+pi/6).",
            "is_closed": False,
            "severity_0_to_1": 0.90,
            "evidence_signals": {
                "p2633_failed_gate_present": "integer_node_pattern_formula_exact_or_repaired" in p2633_failed,
                "p2633_integer_nodes_exact": p2633.get("finite_witness", {}).get("declared_integer_node_audit", {}).get("all_declared_integer_nodes_are_formula_zeros"),
            },
            "next_discharge_test": "Exhaust admissible quotient/reindexing maps from the exact zero lattice to the legacy node list without retuning.",
        },
        {
            "blocker": "inverse_hierarchy_role_transfer",
            "why_it_blocks_full_kernel": "A ToE kernel must say whether distant-octave/Wilson-loop roles survive strict heavy-tailed compression or are rejected.",
            "current_evidence": "P2633 computes strict |K(7)|/|K(1)| below one while legacy amplitude-normalized ratio is above one on the audited grid.",
            "is_closed": False,
            "severity_0_to_1": 0.82,
            "evidence_signals": {
                "p2633_failed_gate_present": "legacy_inverse_hierarchy_numerically_retained" in p2633_failed,
                "strict_k7_over_k1": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("strict"),
                "legacy_k7_over_k1": p2633.get("finite_witness", {}).get("inverse_hierarchy_ratio_abs_k7_over_abs_k1", {}).get("legacy_amplitude_normalized"),
            },
            "next_discharge_test": "Find a sourced distance measure/domain where strict compression preserves the role, or mark the role as modified/rejected.",
        },
        {
            "blocker": "strict_core_selector_source_qw2191",
            "why_it_blocks_full_kernel": "A full strict kernel cannot leave the selector/orientation/source bit external or axiom-only.",
            "current_evidence": "Guardrails and P2635 keep selector/source closure false; QW-2191 remains a strict-core obstruction.",
            "is_closed": False,
            "severity_0_to_1": 1.0,
            "evidence_signals": {
                "p2635_selector_source_gate_failed": "selector_source_closure_obtained" in p2635_failed,
                "q_w_2191_discharged_negative_flag": src["P2635_TOE_NEURAL_UNIVERSE_SIGNATURE"].get("negative_export_flags", {}).get("q_w_2191_discharged"),
            },
            "next_discharge_test": "Export an internal symmetry-breaking/selector source or explicitly mark closure as axiom-augmented non-strict.",
        },
        {
            "blocker": "role_bearing_ltotal_dynamical_source",
            "why_it_blocks_full_kernel": "Stable kernels and analogies must become sourced dynamics before they can enter role-bearing L_total terms.",
            "current_evidence": "P2634/P2635 keep L_total closed; P2570-style roughness source obligations remain unsourced.",
            "is_closed": False,
            "severity_0_to_1": 0.74,
            "evidence_signals": {
                "p2635_negative_ltotal": p2635.get("negative_export_flags", {}).get("role_bearing_ltotal_reenabled"),
                "p2625_status": p2625.get("status"),
            },
            "next_discharge_test": "Derive the variational functional, measure, boundary class, and Euler-Lagrange source from nadsoliton dynamics.",
        },
        {
            "blocker": "blind_empirical_frozen_kernel_confirmation",
            "why_it_blocks_full_kernel": "A ToE-like kernel must survive modern-physics holdouts without sector retuning.",
            "current_evidence": "P2635 proposes CMB/LSS and GW/PTA no-retune tests but does not claim empirical confirmation.",
            "is_closed": False,
            "severity_0_to_1": 0.68,
            "evidence_signals": {
                "p2635_empirical_gate_failed": "empirical_modern_physics_confirmation_obtained" in p2635_failed,
                "modern_physics_test_routes": len(p2635.get("modern_physics_test_matrix", [])),
            },
            "next_discharge_test": "Run preregistered frozen-kernel CMB/LSS or GW/PTA holdout against exponential/spline baselines.",
        },
    ]
    return blockers


def full_kernel_decision(blockers: list[dict[str, Any]], src: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2635 = src["P2635_TOE_NEURAL_UNIVERSE_SIGNATURE"]
    score = p2635.get("weighted_toe_signature_score", {}).get("toe_likeness_score_0_to_1_not_probability", 0.0)
    open_blockers = [row for row in blockers if not row["is_closed"]]
    mean_severity = sum(row["severity_0_to_1"] for row in open_blockers) / len(open_blockers) if open_blockers else 0.0
    gates = {
        "toe_symptoms_present": score >= 0.40,
        "strict_kernel_stability_positive": "strict-kernel stability and reproducibility signature" in p2635.get("weighted_toe_signature_score", {}).get("top_positive_axes", []),
        "all_mandatory_blockers_closed": len(open_blockers) == 0,
        "no_selector_source_obstruction": all(row["blocker"] != "strict_core_selector_source_qw2191" or row["is_closed"] for row in blockers),
        "no_bridge_or_role_transfer_obstruction": all(row["is_closed"] for row in blockers if row["blocker"] in {"legacy_to_strict_bridge_completion", "inverse_hierarchy_role_transfer"}),
        "no_beta_phase_obstruction": all(row["is_closed"] for row in blockers if row["blocker"] in {"positive_beta_or_uv_normalization_source", "phase_frequency_node_gauge_certificate"}),
        "no_empirical_holdout_obstruction": all(row["blocker"] != "blind_empirical_frozen_kernel_confirmation" or row["is_closed"] for row in blockers),
    }
    return {
        "gates": gates,
        "open_blocker_count": len(open_blockers),
        "mandatory_blocker_count": len(blockers),
        "mean_open_blocker_severity_0_to_1": mean_severity,
        "toe_symptom_score_0_to_1_not_probability": score,
        "is_full_kernel_now": all(gates.values()),
        "professorial_verdict": (
            "Strict przejawia silne objawy ToE jako stabilny, jednojądrowy, neural-attention-like kandydat. "
            "Nie jest jednak pełnym kernelem ToE: obecne objawy są koniecznym wsparciem roboczym, ale nie wystarczają bez bridge completion, selector/source, role-transfer, L_total source i blind empirical confirmation."
        ),
        "classification": "ROBUST_TOE_LIKE_WORKING_KERNEL_NOT_FULL_KERNEL",
    }


def minimal_blocker_cut(blockers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sorted_blockers = sorted([row for row in blockers if not row["is_closed"]], key=lambda row: row["severity_0_to_1"], reverse=True)
    return [
        {
            "rank": idx + 1,
            "blocker": row["blocker"],
            "severity_0_to_1": row["severity_0_to_1"],
            "next_discharge_test": row["next_discharge_test"],
        }
        for idx, row in enumerate(sorted_blockers[:4])
    ]


def recommendation(cut: list[dict[str, Any]]) -> str:
    first = cut[0]["blocker"] if cut else "none"
    return (
        "Next honest step: attack the minimal blocker cut rather than adding another global ToE-symptom score. "
        f"Start with `{first}` only if it can be made constructive; operationally the most finite computation remains the phase-frequency/node quotient-exhaustion certificate: enumerate admissible quotient/reindexing maps from the exact cosine zero lattice to the legacy node/gauge list, with no post-hoc retuning. "
        "If no admissible map exists, reject that legacy role and then rerun the role-transfer matrix; if a map exists, use it as the first typed stability-to-role interface theorem before any ToE closure language."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2636/S1586 current ToE blocker lattice and full-kernel decision audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps full-kernel, bridge-completion, beta-source, phase/node, selector, empirical neural-universe, and L_total dynamics content before adding the blocker lattice.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    decision = payload["full_kernel_decision"]
    lines.extend([
        "",
        "## Full-kernel decision",
        "",
        f"ToE symptom score from P2635, not truth probability: `{decision['toe_symptom_score_0_to_1_not_probability']:.4f}`.",
        f"Open mandatory blockers: `{decision['open_blocker_count']}/{decision['mandatory_blocker_count']}`.",
        f"Mean open-blocker severity: `{decision['mean_open_blocker_severity_0_to_1']:.4f}`.",
        f"Is strict a full kernel now? `{decision['is_full_kernel_now']}`.",
        "",
        decision["professorial_verdict"],
        "",
        "## Current blockers",
        "",
        "| blocker | closed? | severity | next discharge test |",
        "| --- | --- | ---: | --- |",
    ])
    for row in payload["current_blocker_lattice"]:
        lines.append(f"| {row['blocker']} | {row['is_closed']} | {row['severity_0_to_1']:.2f} | {row['next_discharge_test']} |")
    lines.extend([
        "",
        "## Minimal blocker cut",
        "",
        "| rank | blocker | severity |",
        "| ---: | --- | ---: |",
    ])
    for row in payload["minimal_blocker_cut"]:
        lines.append(f"| {row['rank']} | {row['blocker']} | {row['severity_0_to_1']:.2f} |")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
        "No full strict kernel, ToE closure, bridge completion, role-transfer, role-bearing `L_total`, or `QW-2191` discharge is claimed.",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    src = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    blockers = current_blocker_lattice(src)
    decision = full_kernel_decision(blockers, src)
    cut = minimal_blocker_cut(blockers)
    payload: dict[str, Any] = {
        "status": "P2636_CURRENT_TOE_BLOCKER_LATTICE_FULL_KERNEL_DECISION_AUDIT_NO_FULL_KERNEL_NO_TOE_CLOSURE",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": source_fingerprints(),
        "current_blocker_lattice": blockers,
        "full_kernel_decision": decision,
        "minimal_blocker_cut": cut,
        "recommended_next_honest_step": recommendation(cut),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)

    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2636/S1586 current ToE blocker lattice and full-kernel decision guard",
        "\n## P2636/S1586 current ToE blocker lattice and full-kernel decision guard\n\n"
        "`P2636/S1586` answers the full-kernel question directly: strict shows real ToE symptoms as a stable, one-kernel, neural-attention-like working architecture, but it is not yet a full ToE kernel.  The current blocker lattice keeps open bridge completion, beta/Z_beta source, phase-frequency node/gauge certificate, inverse-hierarchy role transfer, strict-core selector/QW-2191 source, role-bearing `L_total` dynamics, and blind frozen-kernel empirical confirmation.  Therefore ToE symptoms count as support, not finality.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2636/S1586 full-kernel Ltotal decision guard",
        "\n## P2636/S1586 full-kernel Ltotal decision guard\n\n"
        "`P2636/S1586` blocks the inference from stable/neural strict kernel to full role-bearing dynamics: until the blocker lattice is discharged by typed source theorems and blind empirical tests, the strict kernel remains a robust working kernel rather than a full `L_total`-ready ToE kernel.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "is_full_kernel_now": result["full_kernel_decision"]["is_full_kernel_now"],
        "classification": result["full_kernel_decision"]["classification"],
        "open_blocker_count": result["full_kernel_decision"]["open_blocker_count"],
        "minimal_blocker_cut": result["minimal_blocker_cut"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
