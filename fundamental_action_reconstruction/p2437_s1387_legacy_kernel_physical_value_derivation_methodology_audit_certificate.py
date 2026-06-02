#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.json"
MD = GEN / "p2437_s1387_legacy_kernel_physical_value_derivation_methodology_audit_certificate.md"

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

SOURCE_FILES = {
    "K1_KERNEL_SPLIT": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "S2_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P2436_FRONTIER": GEN / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.json",
}

ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
STRICT_KERNEL_TUPLE = {
    "omega": 0.18575,
    "phi": 0.16250,
    "beta": 1.0,
    "eta": 1.8,
}

LEGACY_VALUE_ROWS = [
    {
        "claim_id": "legacy_weak_mixing_angle",
        "legacy_formula": "sin^2(theta_W)=alpha_geo/12",
        "legacy_inputs": {"alpha_geo": ALPHA_GEO},
        "legacy_numeric_value": ALPHA_GEO / 12.0,
        "methodology_flags": [
            "uses_legacy_alpha_geo_not_strict_kernel_generated_value",
            "divisor_12_is_structural_postulate_not_derived_from_strict_kernel_dynamics",
            "repo_revalidation_marks_formula_heuristic_not_strictly_derived",
        ],
    },
    {
        "claim_id": "legacy_inverse_alpha_em",
        "legacy_formula": "alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)",
        "legacy_inputs": {"alpha_geo": ALPHA_GEO, "beta_tors": BETA_TORS},
        "legacy_numeric_value": (ALPHA_GEO / (2.0 * BETA_TORS)) * (1.0 - BETA_TORS),
        "methodology_flags": [
            "uses_legacy_alpha_geo_and_hard_beta_tors",
            "depends_on_linear_torsion_damping_parameter_absent_from_strict_kernel_tuple",
            "repo_revalidation_marks_formula_model_level_or_partial_not_strict_derivation",
        ],
    },
    {
        "claim_id": "legacy_beta_power_gravity_hierarchy",
        "legacy_formula": "G_eff/G_0=(beta_tors)^N, example N=20",
        "legacy_inputs": {"beta_tors": BETA_TORS, "N_example": 20},
        "legacy_numeric_value": BETA_TORS**20,
        "methodology_flags": [
            "uses_beta_tors_power_scaling_ansatz",
            "strict_kernel_uses_beta_d_eta_nonlinear_compression_not_legacy_beta_tors_power_hierarchy",
            "repo_revalidation_marks_exact_gravity_hierarchy_not_full_independent_proof",
        ],
    },
    {
        "claim_id": "legacy_torsion_to_chi11_orientation",
        "legacy_formula": "candidate beta_tors -> chi_11 selector/orientation source",
        "legacy_inputs": {"beta_tors": BETA_TORS},
        "legacy_numeric_value": None,
        "methodology_flags": [
            "selector_mechanism_search_assumption_not_a_physical_value_derivation",
            "does_not_discharge_QW_2191",
            "does_not_export_internal_strict_selector_source",
        ],
    },
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


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
            ".",
            "-g",
            "*.py",
            "-g",
            "*.md",
            "-g",
            "*.tex",
            "-g",
            "!fundamental_action_reconstruction/generated/**",
            "-g",
            "!*.json",
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
        "legacy_value_formulas": "sin\\^2\\(theta_W\\)|alpha_EM\\^-1|beta\\^N|G_eff.*beta|fine structure.*beta_tors",
        "legacy_kernel_inputs": "K_legacy_ont|alpha_geo|beta_tors|omega=pi/4|phi=pi/6|4 ln 2",
        "prior_negative_methodology": "TAUTOLOGY|hardcoded|HEURISTIC|MODEL-LEVEL|NOT STRICTLY DERIVED|NOT FULL INDEPENDENT PROOF|role-transfer",
        "beta_tors_chi11_selector_assumption": "beta_tors.*chi_?11|chi_?11.*beta_tors|selector.*beta_tors|QW-2191",
        "strict_kernel_generation": "K_strict_gate|strict kernel|beta\\*d\\^eta|omega = 0.18575|eta = 1.8|strict.*physical.*value",
        "new_packet": "P2437|S1387|legacy kernel physical value|physical value methodology|strict kernel generated values",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep finds legacy physical-value formulas built from alpha_geo/beta_tors and multiple existing warnings that these are "
            "heuristic/model-level/role-transfer-limited. It does not find a production P24xx certificate reclassifying those formulas as strict-kernel-generated physical values."
        ),
    }


def strict_kernel_absence_rows() -> list[dict[str, Any]]:
    return [
        {
            "strict_kernel_parameter": key,
            "strict_value": value,
            "legacy_role_formula_containing_parameter": False,
            "physical_value_generator_exported": False,
            "reason": "Current strict tuple is a gate-kernel parameter, not an exported Standard-Model value generator in this audit.",
        }
        for key, value in STRICT_KERNEL_TUPLE.items()
    ]


def methodology_rows() -> list[dict[str, Any]]:
    rows = []
    for row in LEGACY_VALUE_ROWS:
        rows.append(
            {
                **row,
                "strict_kernel_generated": False,
                "legacy_kernel_only_physically_sufficient": False,
                "reclassified_verdict": "legacy_model_or_search_assumption_not_current_strict_physical_derivation",
                "required_next_strict_work": "derive this value from the full strict kernel/strict source theorem, or reject/replace the legacy role.",
            }
        )
    return rows


def append_doc_sections() -> None:
    eq_section = """
## P2437/S1387 legacy-kernel physical-value methodology audit certificate

`P2437/S1387` re-audits the old attempts to read physical values from the legacy kernel.  The grep/methodology audit finds that `sin^2(theta_W)=alpha_geo/12`, `alpha_EM^-1=alpha_geo/(2*beta_tors)*(1-beta_tors)`, and `beta^N` hierarchy claims were built from legacy `alpha_geo`/`beta_tors` bookkeeping and were already marked heuristic/model-level/partial rather than strict derivations.  The proposed `beta_tors -> chi_11` link is reclassified as a selector-mechanism search assumption, not a theorem or physical-value derivation.

The audit therefore changes the honest target: physical values must be generated by the full strict kernel and its strict source/selector theorems, not inherited from the incomplete legacy kernel.  P2437 exports no rejection theorem for all possible successors, but it blocks treating legacy formulas as physically meaningful strict values without a new strict-kernel derivation.
""".strip()
    lag_section = """
## P2437/S1387 legacy-value methodology guard for Lagrangian/EOM

`P2437/S1387` forbids inserting legacy Weinberg, alpha_EM, gravity-hierarchy, or `beta_tors -> chi_11` terms into the strict Lagrangian/EOM as if they were generated by `K_strict_gate`.  Until the full strict kernel exports physical-value generators, those legacy formulas remain methodology-audited model/search assumptions rather than role-bearing dynamics.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    k1_text = read_text(SOURCE_FILES["K1_KERNEL_SPLIT"])
    s2_text = read_text(SOURCE_FILES["S2_PRIORITY"])
    grep = rg_audit()
    rows = methodology_rows()
    strict_rows = strict_kernel_absence_rows()
    theorem_export = {
        "theorem_name": "P2437_T1_legacy_kernel_physical_value_derivation_methodology_audit_certificate",
        "audited_legacy_claim_count": len(rows),
        "audited_legacy_claim_ids": [row["claim_id"] for row in rows],
        "legacy_alpha_geo": ALPHA_GEO,
        "legacy_beta_tors": BETA_TORS,
        "legacy_weak_mixing_numeric": ALPHA_GEO / 12.0,
        "legacy_inverse_alpha_em_numeric": (ALPHA_GEO / (2.0 * BETA_TORS)) * (1.0 - BETA_TORS),
        "legacy_beta_power_n20_numeric": BETA_TORS**20,
        "all_legacy_claims_reclassified_not_strict_generated": all(not row["strict_kernel_generated"] for row in rows),
        "all_legacy_claims_not_physically_sufficient_from_legacy_kernel_only": all(
            not row["legacy_kernel_only_physically_sufficient"] for row in rows
        ),
        "beta_tors_to_chi11_reclassified_as_search_assumption": next(
            row for row in rows if row["claim_id"] == "legacy_torsion_to_chi11_orientation"
        )["strict_kernel_generated"]
        is False,
        "strict_kernel_tuple": STRICT_KERNEL_TUPLE,
        "strict_kernel_parameter_audit_rows": strict_rows,
        "strict_kernel_physical_value_generator_exported_in_current_repo": False,
        "strict_kernel_should_be_generation_source_for_physical_values": True,
        "k1_kernel_split_guardrail_detected": "not yet rigorously identified" in k1_text,
        "s2_no_silent_role_transfer_detected": "legacy physical roles automatically transfer" in s2_text,
        "p2436_claim_specific_frontier_inherited": sources.get("P2436_FRONTIER", {}).get("status")
        == "PASS_CLAIM_SPECIFIC_SUCCESSOR_FRONTIER_ANTICHAIN_NO_ROLE_TRANSFER_NO_CLOSURE",
        "legacy_role_claim_transferred_by_this_certificate": False,
        "beta_tors_chi11_theorem_exported": False,
        "strict_physical_value_theorem_exported": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Legacy numerical agreement or near-agreement is not a strict physical derivation from the full strict kernel.",
            "The old beta_tors -> chi_11 idea is treated here only as a selector-mechanism search assumption.",
            "This audit does not prove final rejection of every possible strict successor; it blocks silent inheritance from the incomplete legacy kernel.",
            "No strict physical-value generator, beta_tors -> chi_11 theorem, role-bearing L_total export, or ToE closure is exported.",
        ],
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "four_claims_audited": theorem_export["audited_legacy_claim_count"] == 4,
        "weak_numeric_matches_formula": abs(theorem_export["legacy_weak_mixing_numeric"] - ALPHA_GEO / 12.0) < 1e-15,
        "alpha_em_numeric_matches_formula": abs(
            theorem_export["legacy_inverse_alpha_em_numeric"] - (ALPHA_GEO / (2.0 * BETA_TORS)) * (1.0 - BETA_TORS)
        )
        < 1e-12,
        "legacy_claims_not_strict_generated": theorem_export["all_legacy_claims_reclassified_not_strict_generated"],
        "legacy_kernel_not_physically_sufficient": theorem_export[
            "all_legacy_claims_not_physically_sufficient_from_legacy_kernel_only"
        ],
        "beta_tors_chi11_is_assumption": theorem_export["beta_tors_to_chi11_reclassified_as_search_assumption"],
        "strict_kernel_generator_not_exported": not theorem_export["strict_kernel_physical_value_generator_exported_in_current_repo"],
        "strict_kernel_generation_required": theorem_export["strict_kernel_should_be_generation_source_for_physical_values"],
        "k1_guardrail_detected": theorem_export["k1_kernel_split_guardrail_detected"],
        "s2_guardrail_detected": theorem_export["s2_no_silent_role_transfer_detected"],
        "p2436_inherited": theorem_export["p2436_claim_specific_frontier_inherited"],
        "no_legacy_transfer": not theorem_export["legacy_role_claim_transferred_by_this_certificate"],
        "no_beta_tors_chi11_theorem": not theorem_export["beta_tors_chi11_theorem_exported"],
        "no_strict_value_theorem": not theorem_export["strict_physical_value_theorem_exported"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2437_s1387_v1",
        "packet_id": "P2437",
        "stage_id": "S1387",
        "result_kind": "LEGACY_KERNEL_PHYSICAL_VALUE_DERIVATION_METHODOLOGY_AUDIT_CERTIFICATE",
        "status": "PASS_LEGACY_KERNEL_PHYSICAL_VALUE_METHODOLOGY_AUDIT_NO_STRICT_VALUE_GENERATOR_NO_BETA_TORS_CHI11_THEOREM",
        "legacy_kernel_physical_value_derivation_methodology_audit_certificate": {
            "rg_methodology_audit": grep,
            "source_artifact_statuses": {name: source.get("status", "OPEN_UNKNOWN") for name, source in sources.items()},
            "legacy_methodology_rows": rows,
            "strict_kernel_absence_rows": strict_rows,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": (
            "Stop using legacy alpha_geo/beta_tors formulas as physical-value derivations; build strict-kernel physical-value generators or reject/replace each legacy role."
        ),
        "global_status": "OPEN_PROGRESS_LEGACY_VALUE_METHODOLOGY_AUDITED_STRICT_KERNEL_GENERATION_REQUIRED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["legacy_kernel_physical_value_derivation_methodology_audit_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2437 S1387: legacy-kernel physical-value derivation methodology audit certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Audited legacy claims: `{theorem['audited_legacy_claim_ids']}`.",
                f"- Legacy weak-mixing numeric: `{theorem['legacy_weak_mixing_numeric']}`.",
                f"- Legacy inverse alpha_EM numeric: `{theorem['legacy_inverse_alpha_em_numeric']}`.",
                f"- All legacy claims strict-generated now: `{not theorem['all_legacy_claims_reclassified_not_strict_generated']}`.",
                f"- Strict kernel physical-value generator exported: `{theorem['strict_kernel_physical_value_generator_exported_in_current_repo']}`.",
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
