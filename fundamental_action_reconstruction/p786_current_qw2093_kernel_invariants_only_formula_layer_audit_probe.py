#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
QW_JSON = REPO / "material_dowodowy" / "korpus_qw_pozostaly" / "raporty_json"

IN_QW2093 = QW_JSON / "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json"
IN_QW2086 = QW_JSON / "report_qw2086_mz_nonanchor_ew_pole_gate.json"
IN_QW2087 = QW_JSON / "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
IN_QW2063 = QW_JSON / "report_qw2063_derivational_reconstruction_shared_flavor_basis.json"
IN_QW2064 = QW_JSON / "report_qw2064_micro_derived_renormalization_constants_gate.json"
IN_ALPHA_GEO = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_F784 = GENERATED / "f784_current_minimal_strict_sm_gr_bridge_registry_packet.json"
IN_P785 = GENERATED / "p785_current_minimal_strict_sm_gr_bridge_registry_boundary_audit_probe.json"

OUT = GENERATED / "p786_current_qw2093_kernel_invariants_only_formula_layer_audit_probe.json"
OUT_SUMMARY = GENERATED / "p786_current_qw2093_kernel_invariants_only_formula_layer_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("-", " ")
        .replace("->", " ")
        .replace("/", " ")
        .replace("_", " ")
        .replace("+", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_QW2093, IN_QW2086, IN_QW2087, IN_QW2063, IN_QW2064, IN_ALPHA_GEO, IN_F784, IN_P785]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P786",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    qw2093 = load_json(IN_QW2093)
    qw2086 = load_json(IN_QW2086)
    qw2087 = load_json(IN_QW2087)
    qw2063 = load_json(IN_QW2063)
    qw2064 = load_json(IN_QW2064)
    alpha_geo_obj = load_json(IN_ALPHA_GEO)
    f784 = load_json(IN_F784)
    p785 = load_json(IN_P785)

    sin2_formula = (qw2093.get("frozen_formulas") or {}).get("sin2_eff_kernel", "")
    alpha_s_mu0_formula = (qw2093.get("frozen_formulas") or {}).get("alpha_s_boundary_mu0", "")
    alpha_s_alpha0_formula = (qw2093.get("frozen_formulas") or {}).get("alpha_s_boundary_alpha0", "")

    sin2_blockers: list[str] = []
    if "alpha_geo" in sin2_formula:
        sin2_blockers.append("alpha_geo_touchpoint_not_yet_reexpressed_as_explicit_strict_invariant_map")
    if contains_all(sin2_formula, ["beta uv", "delta eta"]):
        sin2_blockers.append("micro_radiative_regime_constants_present")
    if contains_all(
        str((qw2086.get("checks") or {}).get("source_sin2_theta_w_eff", "")),
        ["frozen", "alpha geo", "micro radiative shift relation"],
    ):
        sin2_blockers.append("upstream_source_explicitly_labels_frozen_alpha_geo_relation")
    if not (f784.get("targets") or {}).get("sin2_theta_w_eff", {}).get("promotion_blockers"):
        sin2_blockers.append("f784_missing_expected_sin2_promotion_blocker")

    alpha_s_blockers: list[str] = []
    if normalize(alpha_s_mu0_formula) == "m bottom":
        alpha_s_blockers.append("boundary_scale_uses_mass_chain_object_not_kernel_invariant")
    if contains_all(alpha_s_alpha0_formula, ["ln", "m top", "m bottom", "delta eta"]):
        alpha_s_blockers.append("boundary_coupling_uses_mass_ratio_plus_micro_constant_ansatz")
    if contains_all(
        str((qw2087.get("checks") or {}).get("boundary_source", "")),
        ["frozen", "hierarchy log", "boundary ansatz"],
    ):
        alpha_s_blockers.append("upstream_source_explicitly_labels_hierarchy_log_boundary_ansatz")
    if not contains_all(
        json.dumps(qw2063.get("flags") or {}, ensure_ascii=True),
        ["strict first principles foundational constants derived", "false"],
    ):
        alpha_s_blockers.append("q2063_foundational_constant_status_not_detected")
    else:
        alpha_s_blockers.append("mass_chain_still_not_promoted_to_strict_first_principles")

    shared_findings = {
        "alpha_geo_strict_source_present": alpha_geo_obj.get("object") == "alpha_geo_strict_derived_v1",
        "q2064_ci_warning": bool(qw2064.get("ci_warning")),
        "q2064_verdict": qw2064.get("verdict"),
        "p785_boundary_audit_pass": (p785.get("status") == "P785_CURRENT_MINIMAL_BRIDGE_REGISTRY_BOUNDARY_AUDIT_PASS_WITH_BOUNDARIES_EXPLICIT"),
    }

    sin2_ready = len(sin2_blockers) == 0
    alpha_s_ready = len(alpha_s_blockers) == 0
    overall_ready = sin2_ready and alpha_s_ready and not shared_findings["q2064_ci_warning"]

    if overall_ready:
        status = "P786_CURRENT_QW2093_KERNEL_INVARIANTS_ONLY_FORMULA_LAYER_PASS"
    else:
        status = "P786_CURRENT_QW2093_KERNEL_INVARIANTS_ONLY_FORMULA_LAYER_BLOCKED_NONEXPORT"

    artifact = {
        "stage": "P786",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "qw2093": rel(IN_QW2093),
            "qw2086": rel(IN_QW2086),
            "qw2087": rel(IN_QW2087),
            "qw2063": rel(IN_QW2063),
            "qw2064": rel(IN_QW2064),
            "alpha_geo_strict_source": rel(IN_ALPHA_GEO),
            "f784": rel(IN_F784),
            "p785": rel(IN_P785),
        },
        "targets": {
            "sin2_theta_w_eff": {
                "kernel_invariants_only_ready": sin2_ready,
                "frozen_formula": sin2_formula,
                "upstream_source_label": (qw2086.get("checks") or {}).get("source_sin2_theta_w_eff"),
                "blockers": sin2_blockers,
                "current_honest_status": "strict bridge object retained with semantic-transfer blocker" if not sin2_ready else "kernel-invariants-only ready",
            },
            "alpha_s_boundary_mu0_alpha0": {
                "kernel_invariants_only_ready": alpha_s_ready,
                "frozen_formula": {
                    "mu0": alpha_s_mu0_formula,
                    "alpha0": alpha_s_alpha0_formula,
                },
                "upstream_source_label": (qw2087.get("checks") or {}).get("boundary_source"),
                "blockers": alpha_s_blockers,
                "current_honest_status": "strict boundary object retained with frozen-ansatz blocker" if not alpha_s_ready else "kernel-invariants-only ready",
            },
        },
        "shared_findings": shared_findings,
        "current_honest_reading": [
            "The current QW-2093 sin2 formula is still not kernel-invariants-only: alpha_geo does exist as a strict-side source object, but the current QW-2093 relation still leaves it inside a frozen formula lineage rather than an explicit strict invariant map, and the upstream source explicitly labels a frozen alpha_geo relation.",
            "The current QW-2093 alpha_s boundary is still not kernel-invariants-only because it depends on m_top/m_bottom and is explicitly labeled a frozen hierarchy-log boundary ansatz.",
            "QW-2063 still reports strict_first_principles_foundational_constants_derived=false, so the mass-chain layer cannot yet be silently promoted into a first-principles bridge input.",
            "QW-2064 remains useful but still carries a wide-CI warning, so micro constants should not be oversold as a fully stabilized invariants-only closure layer.",
        ],
        "recommended_next_move": [
            "For sin2: either replace alpha_geo with a strict-side invariant observable map, or export a non-transfer theorem stating that the legacy Weinberg-role touchpoint is not part of the strict bridge.",
            "For alpha_s: attempt a dimensionless strict mass-observable ratio boundary in place of the current m_top/m_bottom GeV-level boundary, or demote the hierarchy-log boundary ansatz to non-export support only.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P786",
        "status": status,
        "as_of": AS_OF,
        "sin2_kernel_invariants_only_ready": sin2_ready,
        "alpha_s_kernel_invariants_only_ready": alpha_s_ready,
        "q2064_ci_warning": shared_findings["q2064_ci_warning"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
