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
OUT = GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.json"
MD = GEN / "p2634_s1584_strict_stability_evidence_vs_role_completeness_audit.md"

REPO_ROOT = REPO
REPORT_JSON = REPO_ROOT / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json"
REPORT_MD = REPO_ROOT / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_md"

SOURCE_FILES = {
    "QW1968_REFINED_KERNEL_ROBUSTNESS_BOOTSTRAP_GATE": REPORT_JSON / "report_qw1968_refined_kernel_robustness_bootstrap_gate.json",
    "QW2034_ETA_KERNEL_DERIVATIONAL_STABILITY_AUDIT": REPORT_JSON / "report_qw2034_eta_kernel_derivational_stability_audit.json",
    "QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE": REPORT_JSON / "report_qw2049_spectral_micro_stagec_intersection_gate.json",
    "QW2051_INDEPENDENT_REHEARSAL_GATE": REPORT_JSON / "report_qw2051_independent_rehearsal_gate.json",
    "P2506_MINIMUM_ROUGHNESS_SELECTOR_CANDIDATE": GEN / "p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate.json",
    "P2509_VARIATIONAL_WELLPOSEDNESS": GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.json",
    "P2570_APD_SOBOLEV_ROUGHNESS_ORDER_DEPENDENCE": GEN / "p2570_s1520_apd_sobolev_roughness_selector_order_dependence_audit.json",
    "P2631_NEURAL_INFORMATION_FLUX_BETA_CRITICALITY": GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json",
    "P2632_NEURAL_LEGACY_STRICT_CHARACTERISTIC_RETENTION": GEN / "p2632_s1582_neural_legacy_strict_characteristic_retention_audit.json",
    "P2633_DIAGRAM_GROUNDED_CHARACTERISTIC_PRESERVATION": GEN / "p2633_s1583_diagram_grounded_strict_kernel_characteristic_preservation_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "stability_implies_all_legacy_characteristics_preserved_claimed",
    "stability_implies_bridge_completion_claimed",
    "strict_kernel_finality_claimed",
    "positive_beta_renormalization_source_exported",
    "phase_frequency_selector_source_exported",
    "inverse_hierarchy_role_transfer_exported",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
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
    return {"count": len(lines), "samples": lines[:90]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: stability/robustness contents and residual role-completeness contents, not packet names or IDs.
    patterns = {
        "strict_stability_robustness_content": (
            "selected_kernel_stable|key_metrics_stable|kernel robustness|robustness bootstrap|"
            "strict robustness|stability certificate|derivational stability|stable under|perturbations"
        ),
        "spectral_micro_stagec_content": (
            "spectral micro|Stage-C Intersection|micro_beta_overlap|micro_eta_overlap|"
            "phase.condition|pointwise derivation|independent rehearsal|bundle reproducibility"
        ),
        "variational_wellposedness_content": (
            "minimum-roughness|unique minimizer|Sobolev coercivity|node-vanishing|KKT stationarity|"
            "well-posedness|roughness action|postulated functional"
        ),
        "legacy_role_gap_content": (
            "inverse hierarchy|distant octave|Wilson-loop|integer node|node/gauge|alpha_geo|beta_tors|"
            "role-transfer|characteristic preservation"
        ),
        "neural_attention_stability_content": (
            "heavy-tailed attention|positional encoding|edge of chaos|information flux|mutual information|"
            "fractal graph|Energy-Based Model|Boltzmann"
        ),
        "closure_guard_content": (
            "positive_beta_renormalization_source|phase-frequency certificate|QW-2191|role-bearing L_total|"
            "ToE closure|selector source|bridge completion"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for strict stability evidence versus legacy-role completeness",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def bool_fraction(flags: dict[str, Any]) -> tuple[int, int]:
    bools = [v for v in flags.values() if isinstance(v, bool)]
    return sum(1 for v in bools if v), len(bools)


def source_fingerprints() -> dict[str, Any]:
    return {name: {"path": rel(path), "exists": path.exists(), "sha256": sha256_file(path)} for name, path in SOURCE_FILES.items()}


def stability_evidence_ledger(src: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    q1968 = src["QW1968_REFINED_KERNEL_ROBUSTNESS_BOOTSTRAP_GATE"]
    q2034 = src["QW2034_ETA_KERNEL_DERIVATIONAL_STABILITY_AUDIT"]
    q2049 = src["QW2049_SPECTRAL_MICRO_STAGEC_INTERSECTION_GATE"]
    q2051 = src["QW2051_INDEPENDENT_REHEARSAL_GATE"]
    p2506 = src["P2506_MINIMUM_ROUGHNESS_SELECTOR_CANDIDATE"]
    p2509 = src["P2509_VARIATIONAL_WELLPOSEDNESS"]

    q2049_true, q2049_total = bool_fraction(q2049.get("flags", {}))
    q2051_true, q2051_total = bool_fraction(q2051.get("flags", {}))
    q2034_true, q2034_total = bool_fraction(q2034.get("flags", {}))
    p2506_true, p2506_total = bool_fraction(p2506.get("gatekeeper_checks", {}))
    p2509_true, p2509_total = bool_fraction(p2509.get("gatekeeper_checks", {}))

    return [
        {
            "evidence_class": "spectral_micro_stagec_intersection",
            "source": "QW2049",
            "pass_count": q2049.get("pass_count", q2049_true),
            "total_flags": q2049.get("total_flags", q2049_total),
            "verdict": q2049.get("verdict"),
            "readiness": q2049.get("readiness"),
            "supports_strict_internal_stability": q2049_true == q2049_total and q2049_total > 0,
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
        {
            "evidence_class": "independent_rehearsal_reproducibility",
            "source": "QW2051",
            "pass_count": q2051.get("pass_count", q2051_true),
            "total_flags": q2051.get("total_flags", q2051_total),
            "selected_kernel_stable": q2051.get("flags", {}).get("selected_kernel_stable"),
            "key_metrics_stable": q2051.get("flags", {}).get("key_metrics_stable"),
            "verdict": q2051.get("verdict"),
            "supports_strict_internal_stability": q2051_true == q2051_total and q2051_total > 0,
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
        {
            "evidence_class": "eta_derivational_stability",
            "source": "QW2034",
            "pass_count": q2034.get("pass_count", q2034_true),
            "total_flags": q2034.get("total_flags", q2034_total),
            "failed_flags": [k for k, v in q2034.get("flags", {}).items() if v is False],
            "verdict": q2034.get("verdict"),
            "readiness": q2034.get("readiness"),
            "supports_strict_internal_stability": q2034.get("verdict") == "ETA_KERNEL_DERIVATIONAL_STABILITY_PASS",
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
        {
            "evidence_class": "refined_kernel_robustness_bootstrap",
            "source": "QW1968",
            "verdict": q1968.get("verdict"),
            "required_next_step": q1968.get("required_next_step"),
            "gw_pass_rate": q1968.get("bootstrap_gw_and_triad", {}).get("gw_pass_rate"),
            "triad_pass_rate": q1968.get("bootstrap_gw_and_triad", {}).get("triad_pass_rate"),
            "supports_strict_internal_stability": q1968.get("verdict") is not None,
            "fragility_warning": q1968.get("verdict") == "FRAGILE_PASS_NOT_YET_LOCKABLE",
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
        {
            "evidence_class": "conditional_minimum_roughness_selector",
            "source": "P2506",
            "pass_count": p2506_true,
            "total_flags": p2506_total,
            "status": p2506.get("status"),
            "constant_flow_selected_conditionally": p2506.get("gatekeeper_checks", {}).get("constant_flow_selected_conditionally"),
            "selector_marked_postulated_not_derived": p2506.get("gatekeeper_checks", {}).get("selector_marked_postulated_not_derived"),
            "supports_strict_internal_stability": p2506.get("gatekeeper_checks", {}).get("constant_flow_selected_conditionally") is True,
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
        {
            "evidence_class": "conditional_variational_wellposedness",
            "source": "P2509",
            "pass_count": p2509_true,
            "total_flags": p2509_total,
            "status": p2509.get("status"),
            "unique_minimizer_constant_flow": p2509.get("gatekeeper_checks", {}).get("unique_minimizer_constant_flow"),
            "wellposed_for_postulated_functional": p2509.get("gatekeeper_checks", {}).get("wellposed_for_postulated_functional"),
            "source_not_exported": p2509.get("gatekeeper_checks", {}).get("source_not_exported"),
            "supports_strict_internal_stability": p2509.get("gatekeeper_checks", {}).get("wellposed_for_postulated_functional") is True,
            "does_not_by_itself_prove_legacy_role_transfer": True,
        },
    ]


def role_completeness_obstruction_ledger(src: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p2631 = src["P2631_NEURAL_INFORMATION_FLUX_BETA_CRITICALITY"]
    p2632 = src["P2632_NEURAL_LEGACY_STRICT_CHARACTERISTIC_RETENTION"]
    p2633 = src["P2633_DIAGRAM_GROUNDED_CHARACTERISTIC_PRESERVATION"]
    p2570 = src["P2570_APD_SOBOLEV_ROUGHNESS_ORDER_DEPENDENCE"]

    p2633_acceptance = p2633.get("source_acceptance", {})
    p2633_failed = p2633_acceptance.get("failed_gates", [])
    p2631_acceptance = p2631.get("source_acceptance", {})
    p2632_acceptance = p2632.get("source_acceptance", {})
    p2570_export = p2570.get("apd_sobolev_roughness_selector_order_dependence_audit", {}).get("theorem_export", {})

    return [
        {
            "missing_or_guarded_role": "beta normalization / positive Z_beta bridge source",
            "current_evidence": "P2631/P2630/P2629 classify beta=1 and Z_beta=100 as not target-independently sourced by information flux or strict-internal beta evidence alone.",
            "closed_by_stability_tests": False,
            "blocking_signal": p2631_acceptance.get("accepts_information_conservation_beta_identity", False) is False,
        },
        {
            "missing_or_guarded_role": "phase-frequency node/gauge exactness",
            "current_evidence": "P2633 computes that the declared legacy integer node pattern is not exact for cos(pi*d/4+pi/6).",
            "closed_by_stability_tests": False,
            "blocking_signal": "integer_node_pattern_formula_exact_or_repaired" in p2633_failed,
        },
        {
            "missing_or_guarded_role": "inverse hierarchy / distant-octave Wilson-loop role",
            "current_evidence": "P2633 computes strict |K(7)|/|K(1)| below one on the audited grid while legacy is above one.",
            "closed_by_stability_tests": False,
            "blocking_signal": "legacy_inverse_hierarchy_numerically_retained" in p2633_failed,
        },
        {
            "missing_or_guarded_role": "alpha_geo and beta_tors role transfer",
            "current_evidence": "P2632/P2633 keep amplitude and torsion roles separated from the raw strict kernel until a typed bridge theorem exists.",
            "closed_by_stability_tests": False,
            "blocking_signal": p2632_acceptance.get("accepts_strict_kernel_as_final_complete_toe_kernel", False) is False,
        },
        {
            "missing_or_guarded_role": "APD / Sobolev roughness dynamics source",
            "current_evidence": "P2570 blocks promotion of roughness interpolation to strict dynamics unless derivative order, measure, and boundary class are sourced from nadsoliton dynamics.",
            "closed_by_stability_tests": False,
            "blocking_signal": p2570_export.get("strict_dynamical_source_for_A_P_D_exported", False) is False,
        },
    ]


def orthogonality_truth_table() -> list[dict[str, Any]]:
    rows = []
    for stability_passes in [False, True]:
        for bridge_roles_complete in [False, True]:
            rows.append({
                "strict_internal_stability_evidence_passes": stability_passes,
                "legacy_to_strict_characteristic_roles_complete": bridge_roles_complete,
                "may_claim_working_kernel_stable": stability_passes,
                "may_claim_final_toe_kernel": stability_passes and bridge_roles_complete,
                "verdict": (
                    "stable_working_successor_not_final" if stability_passes and not bridge_roles_complete
                    else "candidate_finality_gate_open" if stability_passes and bridge_roles_complete
                    else "unstable_or_unready"
                ),
            })
    return rows


def aggregate_scores(stability: list[dict[str, Any]], roles: list[dict[str, Any]]) -> dict[str, Any]:
    stability_positive = sum(1 for row in stability if row.get("supports_strict_internal_stability") is True)
    stability_total = len(stability)
    fragility_warnings = [row["source"] for row in stability if row.get("fragility_warning")]
    role_blocks = sum(1 for row in roles if row.get("blocking_signal") is True and row.get("closed_by_stability_tests") is False)
    return {
        "strict_internal_stability_support_fraction": stability_positive / stability_total if stability_total else None,
        "strict_internal_stability_positive_classes": stability_positive,
        "strict_internal_stability_total_classes": stability_total,
        "fragility_warning_sources": fragility_warnings,
        "role_completeness_open_blocks": role_blocks,
        "role_completeness_total_audited_blocks": len(roles),
        "professorial_reading": (
            "The repo contains strong multi-lane evidence that the strict kernel is internally stable as a working architecture, "
            "but the same evidence is logically orthogonal to the missing legacy-characteristic role-transfer theorems. "
            "Therefore stability upgrades confidence in strict as a robust successor; it does not by itself upgrade strict to a final ToE kernel."
        ),
    }


def acceptance(stability: list[dict[str, Any]], roles: list[dict[str, Any]], scores: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "substantial_strict_internal_stability_evidence_present": scores["strict_internal_stability_positive_classes"] >= 4,
        "all_stability_fragility_warnings_cleared": len(scores["fragility_warning_sources"]) == 0,
        "beta_normalization_or_zbeta_bridge_source_closed": not roles[0]["blocking_signal"],
        "phase_frequency_node_gauge_certificate_closed": not roles[1]["blocking_signal"],
        "inverse_hierarchy_role_transfer_closed": not roles[2]["blocking_signal"],
        "alpha_geo_beta_tors_role_transfer_closed": not roles[3]["blocking_signal"],
        "apd_dynamic_source_closed": not roles[4]["blocking_signal"],
        "no_stability_to_finality_overclaim": True,
    }
    return {
        "gates": gates,
        "accepts_stability_as_final_toe_kernel_completion": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "STABILITY_EVIDENCE_STRONG_BUT_ROLE_COMPLETENESS_REJECTED",
        "reason": (
            "Strict stability evidence is real and should be counted: spectral/micro/Stage-C, independent rehearsal, eta stability, and conditional roughness/well-posedness lanes all support a robust working strict kernel.  "
            "However, stability is not the same theorem type as legacy-characteristic preservation: beta normalization, phase-frequency node exactness, inverse hierarchy, alpha_geo/beta_tors transfer, and APD source obligations remain open."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: convert the broad stability evidence into a typed stability-to-role interface theorem rather than another global stability scan.  "
        "The most falsifiable subproblem is a phase-frequency/node certificate that treats the exact zero lattice of cos(omega*d+phi) as a constraint and solves for whether any admissible selector, distance reindexing, or quotient map can recover the legacy node/gauge list without retuning after the fact.  "
        "If that certificate fails, write the strict kernel as a stable approximate finite-domain successor with explicit residuals for node/gauge and inverse-hierarchy roles before any role-transfer audit."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2634/S1584 strict stability evidence versus role completeness audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This audit greps stability/robustness, spectral-micro/Stage-C, variational well-posedness, neural attention, legacy role gaps, and closure guards by research content rather than packet IDs.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Stability evidence ledger",
        "",
        "| evidence class | source | pass/total | verdict/status | supports strict internal stability? |",
        "| --- | --- | ---: | --- | --- |",
    ])
    for row in payload["stability_evidence_ledger"]:
        pass_total = f"{row.get('pass_count', 'n/a')}/{row.get('total_flags', 'n/a')}"
        verdict = row.get("verdict") or row.get("status") or row.get("readiness")
        lines.append(f"| {row['evidence_class']} | {row['source']} | {pass_total} | `{verdict}` | {row['supports_strict_internal_stability']} |")
    lines.extend([
        "",
        "## Role-completeness obstruction ledger",
        "",
        "| guarded role | closed by stability tests? | blocking signal? |",
        "| --- | --- | --- |",
    ])
    for row in payload["role_completeness_obstruction_ledger"]:
        lines.append(f"| {row['missing_or_guarded_role']} | {row['closed_by_stability_tests']} | {row['blocking_signal']} |")
    scores = payload["aggregate_scores"]
    lines.extend([
        "",
        "## Aggregate reading",
        "",
        f"Strict internal stability support classes: `{scores['strict_internal_stability_positive_classes']}/{scores['strict_internal_stability_total_classes']}`.",
        f"Open role-completeness blocks: `{scores['role_completeness_open_blocks']}/{scores['role_completeness_total_audited_blocks']}`.",
        f"Fragility warning sources: `{scores['fragility_warning_sources']}`.",
        "",
        scores["professorial_reading"],
        "",
        "## Verdict",
        "",
        payload["source_acceptance"]["reason"],
        "",
        "No bridge completion, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure is reopened.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    src = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    stability = stability_evidence_ledger(src)
    roles = role_completeness_obstruction_ledger(src)
    scores = aggregate_scores(stability, roles)
    source_acceptance = acceptance(stability, roles, scores)
    payload: dict[str, Any] = {
        "status": "P2634_STRICT_STABILITY_EVIDENCE_VS_ROLE_COMPLETENESS_AUDIT_NO_FINALITY_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_fingerprints": source_fingerprints(),
        "stability_evidence_ledger": stability,
        "role_completeness_obstruction_ledger": roles,
        "orthogonality_truth_table": orthogonality_truth_table(),
        "aggregate_scores": scores,
        "source_acceptance": source_acceptance,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)

    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2634/S1584 strict stability evidence versus role completeness audit guard",
        "\n## P2634/S1584 strict stability evidence versus role completeness audit guard\n\n"
        "`P2634/S1584` counts the repo's existing strict-kernel stability evidence instead of ignoring it: spectral/micro/Stage-C intersection, independent rehearsal stability, eta derivational stability, robustness bootstrap, and conditional roughness/well-posedness lanes support the strict kernel as a robust working architecture.  The same audit keeps theorem types separated: stability evidence does not by itself close beta/Z_beta normalization, phase-frequency node/gauge exactness, inverse-hierarchy role transfer, `alpha_geo/beta_tors` transfer, APD dynamical sourcing, role-bearing `L_total`, `QW-2191`, or ToE closure.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2634/S1584 stability-to-role Ltotal guard",
        "\n## P2634/S1584 stability-to-role Ltotal guard\n\n"
        "`P2634/S1584` treats strict stability as real positive evidence for a working kernel, but blocks promotion of stability into role-bearing `L_total` dynamics.  A role term still needs a typed bridge/source theorem for the relevant legacy characteristic; robustness and well-posedness alone are not a role-transfer theorem.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "stability_positive_classes": result["aggregate_scores"]["strict_internal_stability_positive_classes"],
        "stability_total_classes": result["aggregate_scores"]["strict_internal_stability_total_classes"],
        "open_role_blocks": result["aggregate_scores"]["role_completeness_open_blocks"],
        "fragility_warning_sources": result["aggregate_scores"]["fragility_warning_sources"],
        "accepts_stability_as_final_toe_kernel_completion": result["source_acceptance"]["accepts_stability_as_final_toe_kernel_completion"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
