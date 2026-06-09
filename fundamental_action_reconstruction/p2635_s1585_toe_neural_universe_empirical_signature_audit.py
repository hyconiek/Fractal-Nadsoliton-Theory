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
OUT = GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.json"
MD = GEN / "p2635_s1585_toe_neural_universe_empirical_signature_audit.md"

REPO_ROOT = REPO
REPORT_JSON = REPO_ROOT / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json"
REPORT_MD = REPO_ROOT / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_md"

SOURCE_FILES = {
    "P2634_STABILITY_VS_ROLE_COMPLETENESS": GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json",
    "P2633_DIAGRAM_GROUNDED_RETENTION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
    "P2631_NEURAL_BETA_CRITICALITY": GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json",
    "P2506_MINIMUM_ROUGHNESS_SELECTOR": GEN / "p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate.json",
    "P2509_VARIATIONAL_WELLPOSEDNESS": GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.json",
    "P2570_APD_ROUGHNESS_SOURCE_OBSTRUCTION": GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.json",
    "QW2049_SPECTRAL_MICRO_STAGEC": REPORT_JSON / "report_qw2049_spectral_micro_stagec_intersection_gate.json",
    "QW2051_INDEPENDENT_REHEARSAL": REPORT_JSON / "report_qw2051_independent_rehearsal_gate.json",
    "QW1968_ROBUSTNESS_BOOTSTRAP": REPORT_JSON / "report_qw1968_refined_kernel_robustness_bootstrap_gate.json",
    "QW2076_EMPIRICAL_PREREGISTRATION_MD": REPORT_MD / "RAPORT_QW2076_EMPIRICAL_PREDICTION_PREREGISTRATION.md",
    "FIN_NEURAL_STATUS_MD": REPO_ROOT / "material_dowodowy" / "lean_fin_dowody" / "md" / "FIN_Theory_Status_Report.md",
    "FIN_CLASSICAL_PHYSICS_MD": REPO_ROOT / "material_dowodowy" / "lean_fin_dowody" / "md" / "FIN_VS_FIZYKA_KLASYCZNA.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "toe_closure_claimed",
    "universe_self_learning_proven_as_physical_fact",
    "modern_physics_empirical_confirmation_claimed",
    "strict_kernel_finality_claimed",
    "legacy_role_transfer_revalidated",
    "phase_frequency_selector_source_exported",
    "positive_beta_renormalization_source_exported",
    "q_w_2191_discharged",
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


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


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
    return {"count": len(lines), "samples": lines[:100]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: neural-universe, ToE-signature, and empirical-test contents, not packet numbers.
    patterns = {
        "single_foundation_unification_content": (
            "single kernel|triple sector|common structural function|unification|identity-unification|"
            "single foundation|universal kernel|one fundamental|jednej zasady|pojedyncze.*pole"
        ),
        "neural_universe_self_learning_content": (
            "self-learning|samoucząca|neural universe|sieć neuron|Graph Neural Network|"
            "Energy-Based Model|Boltzmann|attention|positional encoding|gradient descent|learning"
        ),
        "variational_rg_learning_content": (
            "delta S|δS|variational|minimum-roughness|unique minimizer|RG flow|renormalization|"
            "stationarity|Euler-Lagrange|well-posedness"
        ),
        "modern_physics_empirical_content": (
            "CMB|cosmic microwave|large scale structure|LSS|gravitational wave|GW|PTA|LIGO|Virgo|"
            "CKM|PMNS|holdout|preregistration|external confirmatory"
        ),
        "strict_stability_generalization_content": (
            "selected_kernel_stable|key_metrics_stable|Stage-C|spectral micro|derivational stability|"
            "robustness bootstrap|external stress|bundle reproducibility"
        ),
        "nonclosure_obstruction_content": (
            "QW-2191|selector source|positive_beta_renormalization_source|role-transfer|"
            "role-bearing L_total|ToE closure|not supported|inconclusive|fragile"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for ToE signatures, neural-universe symptoms, and empirical-test routes",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def bool_fraction(flags: dict[str, Any]) -> tuple[int, int]:
    bools = [value for value in flags.values() if isinstance(value, bool)]
    return sum(1 for value in bools if value), len(bools)


def source_fingerprints() -> dict[str, Any]:
    return {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()}


def compute_toe_signature_axes(src: dict[str, dict[str, Any]], texts: dict[str, str]) -> list[dict[str, Any]]:
    p2634 = src["P2634_STABILITY_VS_ROLE_COMPLETENESS"]
    p2633 = src["P2633_DIAGRAM_GROUNDED_RETENTION"]
    p2631 = src["P2631_NEURAL_BETA_CRITICALITY"]
    p2506 = src["P2506_MINIMUM_ROUGHNESS_SELECTOR"]
    p2509 = src["P2509_VARIATIONAL_WELLPOSEDNESS"]
    p2570 = src["P2570_APD_ROUGHNESS_SOURCE_OBSTRUCTION"]
    q2049 = src["QW2049_SPECTRAL_MICRO_STAGEC"]
    q2051 = src["QW2051_INDEPENDENT_REHEARSAL"]
    q1968 = src["QW1968_ROBUSTNESS_BOOTSTRAP"]

    p2634_scores = p2634.get("aggregate_scores", {})
    stability_fraction = p2634_scores.get("strict_internal_stability_support_fraction") or 0.0
    role_open = p2634_scores.get("role_completeness_open_blocks", 0)
    role_total = p2634_scores.get("role_completeness_total_audited_blocks", 1) or 1
    role_completion_fraction = max(0.0, 1.0 - role_open / role_total)
    q2049_true, q2049_total = bool_fraction(q2049.get("flags", {}))
    q2051_true, q2051_total = bool_fraction(q2051.get("flags", {}))
    p2506_checks = p2506.get("gatekeeper_checks", {})
    p2509_checks = p2509.get("gatekeeper_checks", {})
    p2570_export = p2570.get("apd_sobolev_roughness_selector_order_dependence_audit", {}).get("theorem_export", {})

    neural_mapping_count = len(p2631.get("neural_prism_mapping", [])) if isinstance(p2631.get("neural_prism_mapping"), list) else 5
    beta_identity_rejected = p2631.get("source_acceptance", {}).get("accepts_information_conservation_beta_identity", False) is False
    p2633_failed = p2633.get("source_acceptance", {}).get("failed_gates", [])
    prereg_text = texts.get("QW2076_EMPIRICAL_PREREGISTRATION_MD", "")
    neural_status_text = texts.get("FIN_NEURAL_STATUS_MD", "") + "\n" + texts.get("FIN_CLASSICAL_PHYSICS_MD", "")

    return [
        {
            "axis": "single-kernel cross-sector unification signature",
            "professorial_reading": "Najmocniejszy objaw ToE: ta sama architektura jądra próbuje obsłużyć kilka sektorów, zamiast doklejać osobne modele sektorowe.",
            "computable_indicators": {
                "qw2049_pass_fraction": q2049_true / q2049_total if q2049_total else None,
                "qw2049_verdict": q2049.get("verdict"),
                "qw2051_selected_kernel_stable": q2051.get("flags", {}).get("selected_kernel_stable"),
            },
            "score_0_to_1_not_probability": 0.85 if q2049_true == q2049_total and q2051.get("flags", {}).get("selected_kernel_stable") else 0.55,
            "weight": 0.18,
            "closure_status": "strong_toe_symptom_not_toe_closure",
        },
        {
            "axis": "strict-kernel stability and reproducibility signature",
            "professorial_reading": "Drugim najmocniejszym objawem ToE jest stabilność: jądro nie rozpada się przy wielu audytach, rehearsal i rekonstrukcjach.",
            "computable_indicators": {
                "p2634_stability_fraction": stability_fraction,
                "fragility_warning_sources": p2634_scores.get("fragility_warning_sources", []),
                "q2051_key_metrics_stable": q2051.get("flags", {}).get("key_metrics_stable"),
            },
            "score_0_to_1_not_probability": max(0.0, 0.90 * stability_fraction - (0.12 if p2634_scores.get("fragility_warning_sources") else 0.0)),
            "weight": 0.16,
            "closure_status": "strong_working_architecture_support_not_role_transfer",
        },
        {
            "axis": "neural architecture signature: positional encoding plus heavy-tailed attention",
            "professorial_reading": "Strict kernel wygląda jak sieciowy operator uwagi: faza pełni rolę kodowania pozycji, a mianownik heavy-tailed attention bias.",
            "computable_indicators": {
                "neural_mapping_fields_count": neural_mapping_count,
                "beta_identity_rejected": beta_identity_rejected,
                "edge_of_chaos_accepts_beta_identity": p2631.get("source_acceptance", {}).get("accepts_edge_of_chaos_beta_identity"),
            },
            "score_0_to_1_not_probability": 0.58 if beta_identity_rejected else 0.75,
            "weight": 0.13,
            "closure_status": "architectural_analogy_positive_but_beta_criticality_unclosed",
        },
        {
            "axis": "self-learning / variational stationarity signature",
            "professorial_reading": "Jeśli Wszechświat jest samouczącą się siecią, formalnym śladem powinien być nie training set, lecz zasada stacjonarności/energia, która wybiera stabilne wagi-jądra.",
            "computable_indicators": {
                "p2506_constant_flow_selected_conditionally": p2506_checks.get("constant_flow_selected_conditionally"),
                "p2506_selector_postulated_not_derived": p2506_checks.get("selector_marked_postulated_not_derived"),
                "p2509_unique_minimizer_constant_flow": p2509_checks.get("unique_minimizer_constant_flow"),
                "p2570_dynamic_source_exported": p2570_export.get("strict_dynamical_source_for_A_P_D_exported"),
            },
            "score_0_to_1_not_probability": 0.48 if p2506_checks.get("constant_flow_selected_conditionally") and p2509_checks.get("unique_minimizer_constant_flow") else 0.25,
            "weight": 0.14,
            "closure_status": "conditional_learning_proxy_source_not_derived",
        },
        {
            "axis": "modern-physics empirical interface signature",
            "professorial_reading": "Teoria zdradza ambicję ToE wtedy, gdy daje mrożone testy dla GW/PTA/flavor/RG/CMB, a nie tylko dopasowuje wewnętrzne artefakty.",
            "computable_indicators": {
                "qw1968_verdict": q1968.get("verdict"),
                "qw1968_gw_pass_rate": q1968.get("bootstrap_gw_and_triad", {}).get("gw_pass_rate"),
                "empirical_preregistration_mentions_holdout": "HOLDOUT" in prereg_text or "holdout" in prereg_text,
            },
            "score_0_to_1_not_probability": 0.42 if q1968.get("verdict") == "FRAGILE_PASS_NOT_YET_LOCKABLE" else 0.55,
            "weight": 0.14,
            "closure_status": "empirical_program_present_but_fragile_or_pending",
        },
        {
            "axis": "legacy-characteristic role completeness signature",
            "professorial_reading": "Największa luka ToE: jądro strict musi jeszcze pokazać, które role legacy naprawdę przetrwały, a które są tylko obrazem roboczym.",
            "computable_indicators": {
                "p2634_role_completion_fraction": role_completion_fraction,
                "p2633_failed_gates": p2633_failed,
            },
            "score_0_to_1_not_probability": max(0.05, 0.25 * role_completion_fraction + 0.05),
            "weight": 0.17,
            "closure_status": "major_open_blocker",
        },
        {
            "axis": "selector/source closure signature",
            "professorial_reading": "Teoria wygląda ToE-podobnie, gdy ma własne selektory, a nie zewnętrzne wybory; tu QW-2191 i źródło fazy/beta nadal blokują finalność.",
            "computable_indicators": {
                "q_w_2191_mentions_present": "QW-2191" in json.dumps(p2634, sort_keys=True),
                "phase_frequency_failed": "integer_node_pattern_formula_exact_or_repaired" in p2633_failed,
                "positive_beta_source_closed": not beta_identity_rejected,
            },
            "score_0_to_1_not_probability": 0.12,
            "weight": 0.10,
            "closure_status": "selector_source_obstruction_open",
        },
        {
            "axis": "neural-universe interpretive continuity signature",
            "professorial_reading": "Repo zawiera język samouczącej się sieci informacyjnej; profesor potraktowałby go jako płodną interpretację, dopóki nie daje ścisłych predykcji niezależnych od fittingu.",
            "computable_indicators": {
                "neural_status_text_present": "neural" in neural_status_text.lower() or "sieć" in neural_status_text.lower(),
                "self_learning_text_present": "samouczą" in neural_status_text.lower() or "self-learning" in neural_status_text.lower(),
            },
            "score_0_to_1_not_probability": 0.50 if neural_status_text else 0.20,
            "weight": 0.08,
            "closure_status": "interpretive_frame_not_empirical_proof",
        },
    ]


def weighted_score(axes: list[dict[str, Any]]) -> dict[str, Any]:
    numerator = sum(axis["score_0_to_1_not_probability"] * axis["weight"] for axis in axes)
    denominator = sum(axis["weight"] for axis in axes)
    score = numerator / denominator if denominator else 0.0
    ranked = sorted(axes, key=lambda row: row["score_0_to_1_not_probability"] * row["weight"], reverse=True)
    blockers = [axis for axis in axes if "open" in axis["closure_status"] or "blocker" in axis["closure_status"] or "unclosed" in axis["closure_status"]]
    return {
        "toe_likeness_score_0_to_1_not_probability": score,
        "top_positive_axes": [row["axis"] for row in ranked[:3]],
        "top_blocking_axes": [row["axis"] for row in blockers],
        "interpretation": (
            "This is a structured symptom score, not a probability that the theory is true.  "
            "The largest ToE signs are single-kernel unification, strict stability/reproducibility, and neural attention architecture; "
            "the largest blockers are selector/source closure and legacy role completeness."
        ),
    }


def modern_physics_test_matrix() -> list[dict[str, Any]]:
    return [
        {
            "test_route": "CMB/LSS fractal-attention correlation test",
            "what_to_measure": "Fit the large-scale two-point correlation tail to a frozen heavy-tailed kernel class and test whether eta near 9/5 is selected without retuning beta/phase.",
            "neural_universe_symptom": "Scale-free memory/attention over cosmic distances rather than a short-range exponential kernel.",
            "pass_condition": "eta and phase constraints survive blind CMB/LSS holdouts better than matched exponential or arbitrary spline baselines.",
            "failure_condition": "The tail exponent or phase has to be retuned per dataset, or exponential/spline baselines dominate.",
        },
        {
            "test_route": "GW/PTA frozen-kernel holdout",
            "what_to_measure": "Use the already frozen strict kernel to predict event-windowed GW/PTA correlation features under detector-coherent controls.",
            "neural_universe_symptom": "The same information-transport kernel generalizes to metric perturbation channels.",
            "pass_condition": "Pre-registered holdout passes with control gap and no kernel retune, clearing the QW1968 fragility warning.",
            "failure_condition": "GW success remains fragile, depends on retuning, or loses against controls.",
        },
        {
            "test_route": "RG / QFT stationarity test",
            "what_to_measure": "Derive the roughness/variational functional from the action, then test whether the constant running-exponent flow is the unique minimizer on independent nodes.",
            "neural_universe_symptom": "Self-learning appears as an energy-based stationary update rule, not as external supervised fitting.",
            "pass_condition": "Derivative order, measure, and boundary class are sourced from nadsoliton dynamics before fitting values.",
            "failure_condition": "Minimum-roughness remains a postulated regularizer only.",
        },
        {
            "test_route": "phase-frequency node/gauge falsification",
            "what_to_measure": "Solve the exact zero lattice of cos(omega*d+phi) and compare it to the claimed gauge/generation node list under explicit quotient/reindexing maps.",
            "neural_universe_symptom": "Discrete tokens/nodes emerge from a positional encoding rather than being manually inserted.",
            "pass_condition": "A single admissible selector or quotient recovers the node list without post-hoc retune.",
            "failure_condition": "The claimed nodes are not zeros and no lawful reindexing exists.",
        },
        {
            "test_route": "cross-sector frozen generalization benchmark",
            "what_to_measure": "Lock the strict kernel once and evaluate flavor, GW/PTA, RG, and cosmology targets against sector-specific baselines.",
            "neural_universe_symptom": "A universal trained representation generalizes across tasks; in physics language, a single action/kernel controls multiple sectors.",
            "pass_condition": "Frozen-kernel performance stays competitive across sectors with uncertainty budgets and adversarial controls.",
            "failure_condition": "Each sector requires its own hidden kernel or normalization convention.",
        },
    ]


def professor_verdict(score: dict[str, Any]) -> str:
    return (
        "Profesorski werdykt: największe oznaki ToE są nie w pojedynczej liczbie beta=1, lecz w zbiegu trzech faktów: "
        "jedna architektura jądra próbuje unifikować sektory, jądro strict jest wielokrotnie stabilne/reprodukowane, "
        "a jego forma ma naturalną interpretację sieciową jako positional encoding plus heavy-tailed attention.  "
        "Jeżeli Wszechświat jest samouczącą się siecią, obecny formalny ślad tej idei jest raczej energy-based/variational: δS=0 lub RG-stacjonarność wybiera wagi-jądra.  "
        "To nadal nie jest domknięcie ToE, bo źródło selektora, beta/Z_beta, node/gauge role i role-transfer legacy->strict pozostają otwarte."
    )


def acceptance(axes: list[dict[str, Any]], score: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "single_kernel_unification_symptom_present": "single-kernel cross-sector unification signature" in score["top_positive_axes"],
        "strict_stability_symptom_present": "strict-kernel stability and reproducibility signature" in score["top_positive_axes"],
        "neural_attention_symptom_present": any("neural" in axis for axis in score["top_positive_axes"]),
        "selector_source_closure_obtained": False,
        "legacy_role_completion_obtained": False,
        "empirical_modern_physics_confirmation_obtained": False,
        "no_toe_overclaim": True,
    }
    return {
        "gates": gates,
        "accepts_toe_closure_or_self_learning_universe_proof": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "TOE_SIGNS_PRESENT_BUT_EMPIRICAL_AND_ROLE_CLOSURE_OPEN",
        "reason": (
            "P2635 finds real ToE-like symptoms: single-kernel unification pressure, robust strict-kernel stability, and a neural/attention interpretation.  "
            "It rejects closure because the modern-physics empirical program is not yet a blind-confirmed source theorem, and the selector/role-transfer obstructions remain theorem-level blockers."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: build the phase-frequency/node certificate as the first stability-to-role interface theorem, then run a frozen-kernel empirical preregistration.  "
        "The immediate computation should solve all admissible quotient/reindexing maps that could reconcile the exact cos(omega*d+phi) zero lattice with the legacy node/gauge list; if none exist, declare the node/gauge role a rejected legacy role.  "
        "In parallel, prepare a no-retune CMB/LSS or GW/PTA holdout protocol where eta=9/5 heavy-tailed attention is compared to exponential and spline baselines."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2635/S1585 ToE neural-universe empirical signature audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps by research content: single-kernel unification, neural/self-learning universe language, variational/RG learning, modern-physics empirical tests, strict stability, and nonclosure obstructions.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Professorial ToE symptom ranking",
        "",
        f"Structured ToE-likeness score, not a truth probability: `{payload['weighted_toe_signature_score']['toe_likeness_score_0_to_1_not_probability']:.4f}`.",
        f"Top positive axes: `{payload['weighted_toe_signature_score']['top_positive_axes']}`.",
        f"Top blocking axes: `{payload['weighted_toe_signature_score']['top_blocking_axes']}`.",
        "",
        payload["professor_verdict"],
        "",
        "| axis | score | closure status |",
        "| --- | ---: | --- |",
    ])
    for axis in payload["toe_signature_axes"]:
        lines.append(f"| {axis['axis']} | {axis['score_0_to_1_not_probability']:.3f} | `{axis['closure_status']}` |")
    lines.extend([
        "",
        "## Neural-universe empirical test matrix",
        "",
        "| route | symptom | pass condition | failure condition |",
        "| --- | --- | --- | --- |",
    ])
    for row in payload["modern_physics_test_matrix"]:
        lines.append(f"| {row['test_route']} | {row['neural_universe_symptom']} | {row['pass_condition']} | {row['failure_condition']} |")
    lines.extend([
        "",
        "## Verdict",
        "",
        payload["source_acceptance"]["reason"],
        "",
        "No ToE closure, self-learning-universe proof, bridge completion, role-transfer, role-bearing `L_total`, or `QW-2191` discharge is claimed.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    src = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    texts = {name: read_text(path) for name, path in SOURCE_FILES.items() if path.suffix != ".json"}
    axes = compute_toe_signature_axes(src, texts)
    score = weighted_score(axes)
    source_acceptance = acceptance(axes, score)
    payload: dict[str, Any] = {
        "status": "P2635_TOE_NEURAL_UNIVERSE_EMPIRICAL_SIGNATURE_AUDIT_NO_TOE_CLOSURE_NO_ROLE_TRANSFER_NO_QW2191",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": source_fingerprints(),
        "toe_signature_axes": axes,
        "weighted_toe_signature_score": score,
        "modern_physics_test_matrix": modern_physics_test_matrix(),
        "professor_verdict": professor_verdict(score),
        "source_acceptance": source_acceptance,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)

    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2635/S1585 ToE neural-universe empirical signature audit guard",
        "\n## P2635/S1585 ToE neural-universe empirical signature audit guard\n\n"
        "`P2635/S1585` gives the professorial ToE/neural-universe reading: the strongest ToE symptoms are single-kernel cross-sector unification pressure, multi-lane strict-kernel stability/reproducibility, and the neural architecture analogy `positional encoding + heavy-tailed attention`.  The self-learning reading is only conditionally visible as an energy-based/variational stationarity proxy (`δS=0`, RG/minimum-roughness), not as a closed physical theorem.  Modern-physics checks must be blind frozen-kernel tests in CMB/LSS, GW/PTA, RG/QFT, and phase-frequency node/gauge channels.  No ToE closure, self-learning-universe proof, role-transfer, role-bearing `L_total`, or `QW-2191` discharge is reopened.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2635/S1585 neural-universe empirical Ltotal guard",
        "\n## P2635/S1585 neural-universe empirical Ltotal guard\n\n"
        "`P2635/S1585` permits the neural-universe analogy only as a test-generating interpretation.  A role-bearing `L_total` term still requires a sourced variational functional and blind empirical interface; attention-language or self-learning language alone is not a dynamical source theorem.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "toe_likeness_score_not_probability": result["weighted_toe_signature_score"]["toe_likeness_score_0_to_1_not_probability"],
        "top_positive_axes": result["weighted_toe_signature_score"]["top_positive_axes"],
        "top_blocking_axes": result["weighted_toe_signature_score"]["top_blocking_axes"],
        "accepts_toe_closure_or_self_learning_universe_proof": result["source_acceptance"]["accepts_toe_closure_or_self_learning_universe_proof"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
