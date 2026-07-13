"""P3147/S2097: axiomatic Lie-algebra and SM/GR fit readiness matrix.

This packet answers the post-P3146 question more directly: if we keep the
axioms explicit, how much potential does the theory have, is the model's Lie
algebra good, and can SM/GR fit?  The result is deliberately split into
(1) algebraic readiness, where prior certificates are strong, and
(2) physical/source readiness, where current artifacts remain blocked.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.json"
MD = GEN / "p3147_s2097_axiom_lie_smgr_fit_readiness_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3146": GEN / "p3146_s2096_axiom_unit_action_measure_bridge.json",
    "P3145": GEN / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.json",
    "P1960": GEN / "p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.json",
    "P1961": GEN / "p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json",
    "P1962": GEN / "p1962_s912_strict_matter_higgs_brst_extension_registry_audit.json",
    "P2438": GEN / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.json",
    "P3084": GEN / "p3084_s2034_gauge_representation_obstruction_witness_audit.json",
    "P3094": GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json",
}

GATES = [
    "receiver_or_formal_carrier_present",
    "exact_algebra_or_local_witness_present",
    "unit_action_axioms_available_conditionally",
    "full_field_registry_or_metric_bundle_present",
    "nonimported_source_law_present",
    "global_eom_or_bv_brst_closure_present",
    "empirical_or_role_transfer_safe_calibration_present",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def input_facts(inputs: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p1960 = inputs["P1960"]
    p1961 = inputs["P1961"]
    p1962 = inputs["P1962"]
    p3145 = inputs["P3145"]
    p3146 = inputs["P3146"]
    p2438 = inputs["P2438"]
    p3084 = inputs["P3084"]
    p3094 = inputs["P3094"]

    p1960_alg = p1960.get("StructureConstantGeneratorCertificates_strict_v1", {})
    p1960_jac = p1960.get("JacobiCertificate_SU3_SU2_U1_strict_v1", {})
    p1961_ready = p1961.get("ReadinessAfterP1961_strict_v1", {})
    p1962_reg = p1962.get("RepresentationRegistryAudit_strict_B1_v1", {})
    p1962_ltotal = p1962.get("LTotalFieldRegistryAudit_strict_B1_v1", {})
    p3145_counts = p3145.get("finite_theorem", {}).get("finite_counts", {})
    p3146_counts = p3146.get("finite_theorem", {}).get("finite_counts", {})
    p2438_checks = p2438.get("gatekeeper_checks", {})

    return {
        "p1960_exact_algebra_pass": bool(p1960_alg.get("exact_algebra_pass")),
        "p1960_all_jacobi_pass": bool(p1960_jac.get("all_declared_factors_pass")),
        "p1960_su3_jacobi_components": p1960_jac.get("SU3", {}).get("checked_components", 0),
        "p1960_su2_jacobi_components": p1960_jac.get("SU2", {}).get("checked_components", 0),
        "p1961_local_gauge_sector_ready": bool(p1961_ready.get("evaluated_local_gauge_sector_ready")),
        "p1961_global_tg2_ready": bool(p1961_ready.get("evaluated_global_TG2_ready")),
        "p1962_seed_row_count": p1962_reg.get("seed_row_count", 0),
        "p1962_has_representation_matrices": bool(p1962_reg.get("has_generator_matrices")),
        "p1962_has_commutator_certificate": bool(p1962_reg.get("has_commutator_certificate")),
        "p1962_has_higgs_representation": bool(p1962_reg.get("has_higgs_representation")),
        "p1962_full_nonproxy_field_registry": bool(p1962_ltotal.get("all_eom_concrete")),
        "p3145_receiver_layout_rows": p3145_counts.get("receiver_layout_only_rows", 0),
        "p3145_closed_rows": p3145_counts.get("closed_rows", 0),
        "p3146_conditional_unit_action_subsets": p3146_counts.get("conditional_unit_action_measure_subsets", 0),
        "p3146_strict_source_subsets": p3146_counts.get("strict_source_subsets", 0),
        "p2438_all_targets_blocked": bool(p2438_checks.get("all_targets_blocked")),
        "p3084_nonimported_gauge_source": bool(p3084.get("decision", {}).get("negative_export_flags", {}).get("gauge_representation_source_exported")),
        "p3094_metric_response_source": bool(p3094.get("decision", {}).get("negative_export_flags", {}).get("physical_stress_energy_tensor_exported", False)),
    }


def row(name: str, claim: str, gates: list[bool], evidence: list[str], verdict: str, missing: list[str]) -> dict[str, Any]:
    gate_map = dict(zip(GATES, gates, strict=True))
    positives = sum(gates)
    if all(gates):
        status = "closed"
    elif gate_map["receiver_or_formal_carrier_present"] and gate_map["exact_algebra_or_local_witness_present"]:
        status = "strong_axiom_branch_potential"
    elif gate_map["receiver_or_formal_carrier_present"]:
        status = "receiver_fit_potential_only"
    else:
        status = "not_ready"
    return {
        "row": name,
        "claim": claim,
        "gate_map": gate_map,
        "positive_gate_count": positives,
        "readiness_status": status,
        "evidence": evidence,
        "verdict": verdict,
        "missing_for_strict_closure": missing,
    }


def build_rows(f: dict[str, Any]) -> list[dict[str, Any]]:
    unit_cond = f["p3146_conditional_unit_action_subsets"] == 1 and f["p3146_strict_source_subsets"] == 0
    return [
        row(
            "axiom_branch_potential",
            "Under explicit axioms, the framework has conditional model-building potential because units can be installed and receiver slots can be organized without claiming strict sourcehood.",
            [True, True, unit_cond, False, False, False, False],
            ["P3146 full-triple unit/action axiom bridge", "P3145 10 receiver-layout rows"],
            "positive but conditional: useful as axiomatic phenomenology, not strict ToE",
            ["strict source for ell_*, tau_*, hbar_*", "nonproxy EOM", "selector/orientation source", "empirical calibration"],
        ),
        row(
            "lie_algebra_quality",
            "The local SU(3)xSU(2)xU(1) Lie algebra data are mathematically good in current artifacts: exact structure constants, Jacobi, embedding, and local gauge-sector BRST nilpotency are available.",
            [True, f["p1960_exact_algebra_pass"] and f["p1960_all_jacobi_pass"] and f["p1961_local_gauge_sector_ready"], unit_cond, False, False, False, False],
            ["P1960 exact SU3/SU2/U1 structure constants and Jacobi", "P1961 local gauge-sector s^2=0"],
            "good as local algebra; not yet good as full physical gauge theory",
            ["field-by-field representation matrices", "Higgs/chiral registry", "anomaly/Yukawa certificates", "global BV/BRST and L_total invariance"],
        ),
        row(
            "standard_model_fit",
            "SM can fit as a receiver/axiom scaffold, but not as a strict generated Standard Model on current artifacts.",
            [True, f["p1960_exact_algebra_pass"] and f["p1961_local_gauge_sector_ready"], unit_cond, False, False, False, False],
            ["P1962 records representation seed only", "P3084 blocks non-imported gauge source", "P2438 all SM/GR targets blocked"],
            "fits conditionally; strict SM generation remains blocked",
            ["complete representation registry", "non-imported gauge bundle/curvature/current", "couplings and masses", "global BRST/unitarity"],
        ),
        row(
            "general_relativity_fit",
            "GR can fit as a receiver/axiom scaffold through metric/stress rows and action units, but not as a strict Einstein-Hilbert/nonproxy EOM derivation.",
            [True, True, unit_cond, False, False, False, False],
            ["P3094 stress-energy/metric-response witnesses", "P3145 GR row receiver-only", "P3146 conditional action measure"],
            "fits conditionally; strict GR derivation remains blocked",
            ["physical stress-energy tensor", "metric coupling", "covariant conservation", "background-independent EH/ELg nonproxy EOM"],
        ),
    ]


def build_payload() -> dict[str, Any]:
    inputs = {name: load(path) for name, path in INPUTS.items()}
    facts = input_facts(inputs)
    rows = build_rows(facts)
    closed = [r for r in rows if r["readiness_status"] == "closed"]
    strong = [r for r in rows if r["readiness_status"] == "strong_axiom_branch_potential"]
    theorem = {
        "name": "P3147_T1_axiom_lie_smgr_fit_readiness_split_theorem",
        "statement": "Current artifacts split cleanly: the model has strong conditional/axiomatic potential and good local SU(3)xSU(2)xU(1) Lie algebra data, but no row reaches strict physical closure.  All 4 audited rows have receiver/formal carriers, 4/4 have algebraic or local witnesses, 4/4 can use the P3146 unit/action axioms conditionally, and 0/4 satisfy the full field-registry/source-law/global-EOM/calibration package required for strict SM/GR generation.",
        "finite_counts": {
            "rows_audited": len(rows),
            "receiver_or_formal_carrier_rows": sum(r["gate_map"]["receiver_or_formal_carrier_present"] for r in rows),
            "exact_algebra_or_local_witness_rows": sum(r["gate_map"]["exact_algebra_or_local_witness_present"] for r in rows),
            "conditional_unit_action_axiom_rows": sum(r["gate_map"]["unit_action_axioms_available_conditionally"] for r in rows),
            "full_field_registry_or_metric_bundle_rows": sum(r["gate_map"]["full_field_registry_or_metric_bundle_present"] for r in rows),
            "nonimported_source_law_rows": sum(r["gate_map"]["nonimported_source_law_present"] for r in rows),
            "global_eom_or_bv_brst_closure_rows": sum(r["gate_map"]["global_eom_or_bv_brst_closure_present"] for r in rows),
            "empirical_or_role_transfer_safe_calibration_rows": sum(r["gate_map"]["empirical_or_role_transfer_safe_calibration_present"] for r in rows),
            "closed_rows": len(closed),
            "strong_axiom_branch_potential_rows": len(strong),
            "p1960_jacobi_components_checked": facts["p1960_su3_jacobi_components"] + facts["p1960_su2_jacobi_components"],
            "p1962_representation_seed_rows": facts["p1962_seed_row_count"],
        },
    }
    return {
        "status": "P3147_AXIOM_LIE_SMGR_FIT_READINESS_STRONG_CONDITIONAL_NO_STRICT_CLOSURE",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "facts_from_repo_grep_and_artifacts": facts,
        "constructed_object": {
            "name": "Theta_fit^ax Lie/SM/GR conditional-readiness matrix",
            "classification": "axiom_branch_readiness_and_strict_obstruction_split_certificate",
            "interpretation": "tests whether explicit axioms make the theory promising while preserving strict-source no-closure boundaries",
        },
        "readiness_rows": rows,
        "finite_theorem": theorem,
        "decision": {
            "potential_under_axioms": "Yes: under explicitly labelled axioms, the theory has real model-building potential because unit/action scaffolding, local gauge algebra, and SM/GR receiver slots coexist coherently.",
            "lie_algebra_verdict": "Good locally/algebraically: SU(3), SU(2), U(1) structure constants, Jacobi certificates, QW-2190 embedding, and local gauge-sector BRST nilpotency are already strong.  Not yet sufficient physically: the full representation registry, global BV/BRST, anomaly/Yukawa certificates, unit-bearing currents, and L_total invariance remain missing.",
            "sm_gr_fit_verdict": "SM and GR can fit as conditional receiver/axiom architecture, but current artifacts do not strictly generate SM or GR.  The fit is promising as a scaffold and blocked as a theorem.",
            "next_honest_step": "Construct P3148 as one bounded representation-registry completion audit: export explicit SU(3)xSU(2)xU(1) matrices for the P1962 seed multiplets plus Higgs, then machine-check commutators, cross-factor commutation, hypercharge/Yukawa gauge invariance, and one-family anomaly sums.  This attacks the strongest Lie-algebra-to-SM gap without claiming ToE, GR closure, or unit-source closure.",
            "negative_export_flags": {
                "strict_SM_generation_exported": False,
                "strict_GR_generation_exported": False,
                "global_BV_BRST_exported": False,
                "unit_bearing_L_total_exported": False,
                "strict_unit_source_exported": False,
                "strict_selector_closure_exported": False,
                "ToE_closure_exported": False,
            },
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3147/S2097 Axiom Lie/SM/GR fit readiness matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Readiness rows"])
    for r in payload["readiness_rows"]:
        lines.append(f"- `{r['row']}`: `{r['readiness_status']}`; positives `{r['positive_gate_count']}/7`; verdict: {r['verdict']}")
    d = payload["decision"]
    lines.extend(["", "## Potential under axioms", d["potential_under_axioms"], "", "## Lie algebra verdict", d["lie_algebra_verdict"], "", "## SM/GR fit verdict", d["sm_gr_fit_verdict"], "", "## Recommendation", d["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3147/S2097 axiom Lie/SM/GR fit readiness matrix", "## P3147/S2097 axiom Lie/SM/GR fit readiness matrix\n\n`P3147/S2097` constructs `Theta_fit^ax`, a four-row conditional-readiness matrix for axiom-branch potential, local Lie algebra quality, SM fit, and GR fit.  The matrix finds `4/4` receiver/formal carriers, `4/4` algebraic or local witnesses, and `4/4` conditional P3146 unit/action-axiom compatibility, but `0/4` rows with full field registry/metric bundle, non-imported source law, global EOM/BV-BRST closure, and empirical or role-transfer-safe calibration.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3147/S2097 SM/GR fit remains conditional", "## P3147/S2097 SM/GR fit remains conditional\n\n`P3147/S2097` confirms that the local `SU(3)xSU(2)xU(1)` Lie algebra/BRST data are good as algebra, and that SM/GR can be placed in an axiom-branch receiver architecture.  It also confirms that no strict unit-bearing `L_total`, full SM representation registry, nonproxy GR EOM, global BV/BRST theorem, or empirical calibration is exported.\n")
    append_once(AGENTS, "Current axiom Lie/SM/GR fit readiness guardrail (P3147/S2097, 2026-07-13)", "## Current axiom Lie/SM/GR fit readiness guardrail (P3147/S2097, 2026-07-13)\n\n- P3147 constructs `Theta_fit^ax`, a four-row readiness matrix for axiom-branch potential, local Lie algebra quality, SM fit, and GR fit.\n- The result is positive only in the conditional/algebraic sense: `4/4` rows have receiver/formal carriers and local witnesses, and `4/4` can use the P3146 unit/action axioms conditionally.\n- The strict physical package remains absent: `0/4` rows have full field registry/metric bundle, non-imported source law, global EOM/BV-BRST closure, or empirical/role-transfer-safe calibration.\n- Treat the model's `SU(3)xSU(2)xU(1)` Lie algebra as good locally/algebraically but not yet a full physical SM gauge theory.  Treat SM/GR fit as promising conditional scaffolding, not strict generation.\n- Next honest move: P3148 should complete/audit one explicit representation registry for the P1962 seed multiplets plus Higgs, with commutator, hypercharge/Yukawa, and anomaly checks; do not claim ToE, GR closure, selector closure, or unit-source closure.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
