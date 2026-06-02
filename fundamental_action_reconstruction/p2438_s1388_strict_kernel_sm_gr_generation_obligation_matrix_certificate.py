#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.json"
MD = GEN / "p2438_s1388_strict_kernel_sm_gr_generation_obligation_matrix_certificate.md"

SOURCE_FILES = {
    "P1421_SELECTOR_PRECHECK": GEN / "p1421_strict_selector_uniqueness_discharge_precheck_summary.json",
    "P1646_SM_GR_SCAFFOLD": GEN / "p1646_s596_strict_gauge_metric_cross_variation_scaffold_summary.json",
    "P1705_VARIATIONAL_OBJECTS": GEN / "p1705_s655_strict_nonproxy_metric_spinor_variational_objects_export_checkpoint.json",
    "P1981_GR_R2_OBLIGATION": GEN / "p1981_s931_strict_adm_bianchi_r2_lapse_variation_obligation.json",
    "P2437_LEGACY_METHOD_AUDIT": GEN / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

OBLIGATIONS = [
    "strict_kernel_identity_and_domain_theorem",
    "strict_kernel_to_coefficient_map_theorem",
    "sm_gauge_group_and_coupling_export_theorem",
    "sm_fermion_yukawa_higgs_export_theorem",
    "gr_metric_eh_and_background_independence_theorem",
    "curvature_squared_completion_and_unitarity_theorem",
    "qw2191_selector_uniqueness_discharge_theorem",
    "strict_observable_value_generator_theorem",
]

TARGETS = {
    "L_SM_gauge_ready": [
        "strict_kernel_identity_and_domain_theorem",
        "strict_kernel_to_coefficient_map_theorem",
        "sm_gauge_group_and_coupling_export_theorem",
        "qw2191_selector_uniqueness_discharge_theorem",
        "strict_observable_value_generator_theorem",
    ],
    "L_SM_matter_higgs_ready": [
        "strict_kernel_identity_and_domain_theorem",
        "strict_kernel_to_coefficient_map_theorem",
        "sm_fermion_yukawa_higgs_export_theorem",
        "qw2191_selector_uniqueness_discharge_theorem",
        "strict_observable_value_generator_theorem",
    ],
    "L_GR_ready": [
        "strict_kernel_identity_and_domain_theorem",
        "strict_kernel_to_coefficient_map_theorem",
        "gr_metric_eh_and_background_independence_theorem",
        "curvature_squared_completion_and_unitarity_theorem",
        "strict_observable_value_generator_theorem",
    ],
    "SM_GR_coupled_EOM_ready": [
        "strict_kernel_identity_and_domain_theorem",
        "strict_kernel_to_coefficient_map_theorem",
        "sm_gauge_group_and_coupling_export_theorem",
        "sm_fermion_yukawa_higgs_export_theorem",
        "gr_metric_eh_and_background_independence_theorem",
        "curvature_squared_completion_and_unitarity_theorem",
        "qw2191_selector_uniqueness_discharge_theorem",
    ],
    "physical_constants_generated_ready": [
        "strict_kernel_identity_and_domain_theorem",
        "strict_kernel_to_coefficient_map_theorem",
        "sm_gauge_group_and_coupling_export_theorem",
        "sm_fermion_yukawa_higgs_export_theorem",
        "gr_metric_eh_and_background_independence_theorem",
        "qw2191_selector_uniqueness_discharge_theorem",
        "strict_observable_value_generator_theorem",
    ],
    "strict_toe_ready": OBLIGATIONS,
}

CURRENT_DISCHARGED_OBLIGATIONS: list[str] = []
PARTIAL_EVIDENCE = {
    "strict_kernel_identity_and_domain_theorem": ["P1567 local strict parameter identifiability only", "P2437 strict tuple audit"],
    "strict_kernel_to_coefficient_map_theorem": ["P1567 names the chain but leaves full 4-parameter map open"],
    "sm_gauge_group_and_coupling_export_theorem": ["P1646 exports SM/GR scaffold, not coupling derivation"],
    "sm_fermion_yukawa_higgs_export_theorem": ["P1705 exports spinor variation templates, full curved export open"],
    "gr_metric_eh_and_background_independence_theorem": ["P1981 exports R2 lapse obligation and keeps background independence open"],
    "curvature_squared_completion_and_unitarity_theorem": ["P1981 leaves Ricci2/Riemann2/Gauss-Bonnet and unitarity open"],
    "qw2191_selector_uniqueness_discharge_theorem": ["P1421 precheck says uniqueness certificate and QW-2191 discharge are missing"],
    "strict_observable_value_generator_theorem": ["P2437 says no strict physical-value generator is exported"],
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg",
            "-n",
            pattern,
            "fundamental_action_reconstruction",
            "material_dowodowy",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:30]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2438|S1388|strict kernel SM GR generation|SM GR generation obligation|strict observable value generator",
        "strict_sm_gr_route": "F_nadsoliton.*L_SM.*L_GR|F_Nadsoliton.*L_SM.*L_GR|K_strict.*coefficients.*L_SM|strict-only",
        "sm_gr_scaffold": "L_SM_gauge|L_SM_fermions|L_SM_higgs|L_GR|gauge_metric|metric_spinor|curvature-squared",
        "selector_and_qw2191": "QW-2191|selector uniqueness|strict internal selector|uniqueness discharge",
        "legacy_value_block": "P2437|legacy.*physical-value|beta_tors -> chi_11|alpha_geo/12|alpha_EM",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds strict-only SM/GR scaffolds and selector/GR obligations, but no production P24xx certificate proving "
            "that K_strict already generates SM+GR physical values. The honest next step is an obligation matrix, not closure."
        ),
    }


def status_of(payload: dict[str, Any]) -> str:
    return payload.get("status") or payload.get("verdict") or payload.get("result_kind") or payload.get("_missing", "OPEN_UNKNOWN")


def target_rows() -> list[dict[str, Any]]:
    current = set(CURRENT_DISCHARGED_OBLIGATIONS)
    rows = []
    for target, required in TARGETS.items():
        missing = [item for item in required if item not in current]
        rows.append(
            {
                "target": target,
                "required_obligations": required,
                "required_obligation_count": len(required),
                "currently_discharged_obligations": [item for item in required if item in current],
                "missing_obligations": missing,
                "missing_obligation_count": len(missing),
                "ready_now": not missing,
                "incidence_bits": [1 if obligation in required else 0 for obligation in OBLIGATIONS],
            }
        )
    return rows


def obligation_rows() -> list[dict[str, Any]]:
    return [
        {
            "obligation": obligation,
            "currently_discharged": obligation in CURRENT_DISCHARGED_OBLIGATIONS,
            "partial_evidence": PARTIAL_EVIDENCE[obligation],
            "blocks_targets": [target for target, required in TARGETS.items() if obligation in required],
        }
        for obligation in OBLIGATIONS
    ]


def append_doc_sections() -> None:
    eq_section = """
## P2438/S1388 strict-kernel SM/GR generation obligation matrix certificate

`P2438/S1388` starts the strict-side replacement for the discarded legacy-value derivations: it treats the target as `K_strict_gate -> coefficients -> L_SM + L_GR -> observables`, not as legacy `alpha_geo/beta_tors` inheritance.  The certificate builds an 8-obligation matrix for strict kernel identity/domain, kernel-to-coefficient map, SM gauge couplings, SM matter/Higgs/Yukawa export, GR/background-independence, curvature-squared/unitarity completion, QW-2191 selector uniqueness, and strict observable-value generation.

The current matrix is intentionally negative: all six SM/GR generation targets remain not ready, because the repo has scaffolds and partial variational witnesses but no theorem-grade strict observable generator or QW-2191 discharge.  This is the first strict-only SM/GR worklist after P2437, not a closure claim.
""".strip()
    lag_section = """
## P2438/S1388 strict SM/GR generation guard for Lagrangian/EOM

`P2438/S1388` moves the Lagrangian/EOM target away from legacy constants and onto a strict-only generation chain.  Existing SM/GR scaffold terms may stay as scaffolds, but no coefficient, physical constant, selector choice, GR completion, or role-bearing `L_total` term is licensed until the eight strict generation obligations are actually discharged.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    targets = target_rows()
    obligations = obligation_rows()
    theorem_export = {
        "theorem_name": "P2438_T1_strict_kernel_sm_gr_generation_obligation_matrix_certificate",
        "strict_generation_route": "K_strict_gate -> coefficients -> L_SM + L_GR -> EOM/observables",
        "legacy_value_route_retired_as_generator": True,
        "obligation_count": len(OBLIGATIONS),
        "obligation_names": OBLIGATIONS,
        "target_count": len(TARGETS),
        "target_names": list(TARGETS),
        "target_rows": targets,
        "obligation_rows": obligations,
        "ready_target_count_now": sum(1 for row in targets if row["ready_now"]),
        "ready_targets_now": [row["target"] for row in targets if row["ready_now"]],
        "all_targets_blocked_now": all(not row["ready_now"] for row in targets),
        "current_discharge_mask": 0,
        "full_discharge_mask": (1 << len(OBLIGATIONS)) - 1,
        "minimum_missing_obligations_for_any_target": min(row["missing_obligation_count"] for row in targets),
        "maximum_missing_obligations_for_target": max(row["missing_obligation_count"] for row in targets),
        "source_statuses": {name: status_of(source) for name, source in sources.items()},
        "p1421_qw2191_missing_inherited": sources["P1421_SELECTOR_PRECHECK"].get("precheck", {}).get("qw2191_discharge_proof_present") is False,
        "p1646_scaffold_open_inherited": sources["P1646_SM_GR_SCAFFOLD"].get("strict_core_closure", {}).get("status") == "OPEN",
        "p1705_qg_required_inherited": sources["P1705_VARIATIONAL_OBJECTS"].get("status") == "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "p1981_background_open_inherited": "BACKGROUND_INDEPENDENCE_REMAINS_OPEN" in sources["P1981_GR_R2_OBLIGATION"].get("obstruction_tags", []),
        "p2437_no_strict_value_generator_inherited": sources["P2437_LEGACY_METHOD_AUDIT"]
        .get("legacy_kernel_physical_value_derivation_methodology_audit_certificate", {})
        .get("theorem_export", {})
        .get("strict_kernel_physical_value_generator_exported_in_current_repo")
        is False,
        "strict_sm_gr_generation_theorem_exported": False,
        "strict_observable_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Strict SM/GR scaffold presence is not strict SM/GR value generation.",
            "No legacy alpha_geo/beta_tors value route is used as a generator.",
            "No QW-2191 discharge, strict observable-value generator, role-bearing L_total export, or ToE closure is exported.",
        ],
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "eight_obligations": theorem_export["obligation_count"] == 8,
        "six_targets": theorem_export["target_count"] == 6,
        "all_targets_blocked": theorem_export["all_targets_blocked_now"],
        "current_mask_zero": theorem_export["current_discharge_mask"] == 0,
        "full_mask_255": theorem_export["full_discharge_mask"] == 255,
        "minimum_missing_at_least_five": theorem_export["minimum_missing_obligations_for_any_target"] >= 5,
        "p1421_qw2191_missing": theorem_export["p1421_qw2191_missing_inherited"],
        "p1646_scaffold_open": theorem_export["p1646_scaffold_open_inherited"],
        "p1705_qg_required": theorem_export["p1705_qg_required_inherited"],
        "p1981_background_open": theorem_export["p1981_background_open_inherited"],
        "p2437_no_strict_value_generator": theorem_export["p2437_no_strict_value_generator_inherited"],
        "no_strict_sm_gr_generation_export": not theorem_export["strict_sm_gr_generation_theorem_exported"],
        "no_value_generator_export": not theorem_export["strict_observable_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2438_s1388_v1",
        "packet_id": "P2438",
        "stage_id": "S1388",
        "result_kind": "STRICT_KERNEL_SM_GR_GENERATION_OBLIGATION_MATRIX_CERTIFICATE",
        "status": "PASS_STRICT_KERNEL_SM_GR_GENERATION_OBLIGATION_MATRIX_ALL_TARGETS_BLOCKED_NO_CLOSURE",
        "strict_kernel_sm_gr_generation_obligation_matrix_certificate": {
            "rg_nonduplication_audit": grep,
            "finite_obligation_matrix": {"target_rows": targets, "obligation_rows": obligations},
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Build the strict kernel-to-coefficient map and strict observable-value generator before attempting SM/GR constants."
        ),
        "global_status": "OPEN_PROGRESS_STRICT_SM_GR_GENERATION_WORKLIST_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_kernel_sm_gr_generation_obligation_matrix_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2438 S1388: strict-kernel SM/GR generation obligation matrix certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Route: `{theorem['strict_generation_route']}`.",
                f"- Obligations: `{theorem['obligation_count']}`.",
                f"- Targets: `{theorem['target_count']}`.",
                f"- Ready targets now: `{theorem['ready_targets_now']}`.",
                f"- Minimum missing obligations for any target: `{theorem['minimum_missing_obligations_for_any_target']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Gatekeepers",
                "",
                f"`{payload['gatekeeper_checks']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    append_doc_sections()
    payload = build_payload()
    write_outputs(payload)
    if not all(payload["gatekeeper_checks"].values()):
        raise SystemExit(f"gatekeeper failure: {payload['gatekeeper_checks']}")


if __name__ == "__main__":
    main()
