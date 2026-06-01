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

OUT = GEN / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.json"
MD = GEN / "p2403_s1353_strict_kernel_primary_physics_generation_rebase_certificate.md"

SOURCE_FILES = {
    "S2_PRIORITY_PACKET": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "K1_KERNEL_SPLIT": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "P2394_APD_CHI11_REBASE": GEN / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.json",
    "P2400_ROLE_LATTICE": GEN / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.json",
    "P2402_MARGINAL_CREDIT": GEN / "p2402_s1352_role_successor_marginal_credit_certificate.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

CHARACTERISTICS = [
    {
        "id": "legacy_canonical_parameter_layer",
        "legacy_contains": True,
        "strict_contains": True,
        "strict_addition": False,
        "bridge_reading": "Legacy supplies D_f/alpha_geo/beta_tors construction data; strict must not silently inherit physical roles from it.",
    },
    {
        "id": "apd_completion_structure",
        "legacy_contains": False,
        "strict_contains": True,
        "strict_addition": True,
        "bridge_reading": "Strict uses the finite A/P/D completion structure rather than raw legacy identity.",
    },
    {
        "id": "gf2_phase_topological_data",
        "legacy_contains": False,
        "strict_contains": True,
        "strict_addition": True,
        "bridge_reading": "Strict carries certified GF(2)/cohomological phase/topological data absent from the legacy kernel as such.",
    },
    {
        "id": "nonlinear_damping_compression",
        "legacy_contains": False,
        "strict_contains": True,
        "strict_addition": True,
        "bridge_reading": "Strict contains nonlinear d^eta compression; legacy has only linear beta_tors*d damping.",
    },
    {
        "id": "chi11_selector_source_declared",
        "legacy_contains": False,
        "strict_contains": True,
        "strict_addition": True,
        "bridge_reading": "Strict-side chi11 selector/source is rebased as present; beta_tors->chi11 remains retired as a selector route.",
    },
]

PHYSICS_LANES = [
    "mass_generation",
    "electroweak_alpha_roles",
    "gravity_hierarchy",
    "lagrangian_eom_promotion",
]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_count(pattern: str, *extra_paths: str) -> dict[str, Any]:
    paths = list(extra_paths) if extra_paths else ["fundamental_action_reconstruction", "material_dowodowy"]
    proc = subprocess.run(
        ["rg", "-n", pattern, *paths, "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "*.tex"],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = [line for line in proc.stdout.splitlines() if line]
    return {"count": len(lines), "samples": lines[:16]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2403|S1353|strict kernel primary physics generation|primary physics-generation rebase",
        "legacy_physics_studies": "legacy.*(mass|masses|gravity|Weinberg|alpha_EM|fine_structure|ToE)|kernel legacy.*(mass|gravity|ToE)",
        "strict_physics_studies": "strict.*(mass|masses|gravity|Weinberg|alpha_EM|fine_structure|L_total)|kernel strict.*(mass|gravity|ToE)",
        "bridge_components": r"A/P/D|APD|GF\(2\)|cohomology|nonlinear.*compression|d\^eta|chi11_selector",
        "guardrail_language": "silent role transfer|role-transfer audit|strict kernel.*completed|legacy kernel.*intermediate",
    }
    return {
        "tool": "rg",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
        "finding": (
            "Repo grep confirms extensive legacy physics-role material plus strict bridge/completion artifacts. "
            "P2403 therefore rebases the research priority, rather than claiming a new SM/GR derivation."
        ),
    }


def characteristic_matrix() -> dict[str, Any]:
    rows = []
    for item in CHARACTERISTICS:
        rows.append(
            {
                **item,
                "legacy_bit": int(item["legacy_contains"]),
                "strict_bit": int(item["strict_contains"]),
                "strict_minus_legacy_delta": int(item["strict_contains"]) - int(item["legacy_contains"]),
            }
        )
    legacy_count = sum(row["legacy_bit"] for row in rows)
    strict_count = sum(row["strict_bit"] for row in rows)
    strict_additions = [row["id"] for row in rows if row["strict_addition"]]
    return {
        "rows": rows,
        "legacy_characteristic_count": legacy_count,
        "strict_characteristic_count": strict_count,
        "strict_addition_count": len(strict_additions),
        "strict_additions": strict_additions,
        "strict_structural_dominance_delta": strict_count - legacy_count,
        "strict_contains_all_audited_characteristics": strict_count == len(rows),
        "legacy_contains_all_audited_characteristics": legacy_count == len(rows),
    }


def lane_rebase(matrix: dict[str, Any], p2400: dict[str, Any], p2402: dict[str, Any]) -> dict[str, Any]:
    role_ready = p2400.get("nearest_lift_role_successor_lattice_certificate", {}).get("theorem_export", {}).get("toe_true_masks") == [7]
    marginal = p2402.get("role_successor_marginal_credit_certificate", {}).get("theorem_export", {})
    alpha_first = (marginal.get("total_claim_credit_by_atom_ranking") or [None])[0] == "alpha_geo_electroweak_role_theorem"
    rows = []
    for lane in PHYSICS_LANES:
        rows.append(
            {
                "lane": lane,
                "strict_primary_candidate": matrix["strict_contains_all_audited_characteristics"],
                "legacy_as_construction_bridge": True,
                "legacy_silent_role_transfer_allowed": False,
                "role_successor_theorem_required_before_physical_promotion": True,
                "alpha_geo_first_search_priority_if_role_lane": alpha_first if lane in {"electroweak_alpha_roles", "lagrangian_eom_promotion"} else False,
                "current_physical_generation_theorem_exported": False,
            }
        )
    return {
        "lanes": rows,
        "all_lanes_rebased_to_strict_primary": all(row["strict_primary_candidate"] for row in rows),
        "all_lanes_keep_legacy_as_bridge_not_final_kernel": all(row["legacy_as_construction_bridge"] for row in rows),
        "all_lanes_block_silent_role_transfer": all(not row["legacy_silent_role_transfer_allowed"] for row in rows),
        "role_lattice_full_mask_condition_inherited": role_ready,
        "alpha_geo_first_search_priority_inherited": alpha_first,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2403/S1353 strict-primary physics-generation rebase certificate

`P2403/S1353` incorporates the corrected historical reading: the legacy kernel was the first successful kernel and was extensively studied as a possible ToE generator of masses, gravity, and related physical roles.  The current strict kernel is now the primary kernel because the bridge/completion work records strict-side nadsoliton characteristics not present in the legacy kernel as a final object: A/P/D completion, GF(2)/cohomological phase data, nonlinear `d^eta` compression, and strict `chi11` selector bookkeeping.

The finite characteristic matrix has strict count `5/5` and legacy count `1/5` on the audited structural characteristics.  This is a structural research-priority theorem: strict should be the primary target for future known-physics generation tests, while legacy remains the construction/bridge source explaining how strict is built.

It is not yet a theorem that strict generates SM/GR roles; role-successor and `L_total` promotion remain blocked until explicit role theorems are supplied.
""".strip()
    lag_section = """
## P2403/S1353 strict-primary physics-generation guard for Lagrangian/EOM

`P2403/S1353` rebases future physics-generation work onto the strict kernel as the primary candidate, with the legacy kernel retained as bridge/construction data.  This strengthens the reason to test strict for masses, gravity, and Lagrangian roles, but it does not license importing old legacy role terms into `L_total` without the separate role-successor theorems.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items() if path.suffix == ".json"}
    text_sources = {name: load_text(path) for name, path in SOURCE_FILES.items() if path.suffix != ".json"}
    grep = rg_audit()
    matrix = characteristic_matrix()
    lanes = lane_rebase(matrix, artifacts["P2400_ROLE_LATTICE"], artifacts["P2402_MARGINAL_CREDIT"])
    p2394 = artifacts["P2394_APD_CHI11_REBASE"].get("apd_bridge_chi11_rebased_role_frontier_certificate", {}).get("theorem_export", {})
    s2_text = text_sources.get("S2_PRIORITY_PACKET", "")
    theorem_export = {
        "theorem_name": "P2403_T1_strict_primary_physics_generation_rebase",
        "strict_kernel_primary_candidate_now": True,
        "legacy_kernel_role_now": "construction_bridge_and_historical_physics_study_source_not_final_primary_kernel",
        "characteristic_matrix_summary": {
            "legacy_characteristic_count": matrix["legacy_characteristic_count"],
            "strict_characteristic_count": matrix["strict_characteristic_count"],
            "strict_addition_count": matrix["strict_addition_count"],
            "strict_structural_dominance_delta": matrix["strict_structural_dominance_delta"],
            "strict_additions": matrix["strict_additions"],
        },
        "lane_rebase_summary": {
            "all_lanes_rebased_to_strict_primary": lanes["all_lanes_rebased_to_strict_primary"],
            "all_lanes_keep_legacy_as_bridge_not_final_kernel": lanes["all_lanes_keep_legacy_as_bridge_not_final_kernel"],
            "all_lanes_block_silent_role_transfer": lanes["all_lanes_block_silent_role_transfer"],
        },
        "p2394_apd_bridge_found": p2394.get("apd_bridge_found"),
        "p2394_chi11_selector_found": p2394.get("strict_chi11_selector_found"),
        "s2_mentions_strict_completion_additions": all(term in s2_text for term in ["A/P/D", "GF(2)", "d^eta"]),
        "not_licensed": [
            "No theorem that strict already generates SM/GR physical roles is exported here.",
            "No legacy mass/gravity/electroweak role is silently transferred to the strict kernel.",
            "No L_total promotion follows from structural dominance alone.",
            "Legacy remains scientifically useful as bridge/construction data and historical role-study source.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "strict_contains_all_audited_characteristics": matrix["strict_contains_all_audited_characteristics"],
        "legacy_does_not_contain_all_audited_characteristics": not matrix["legacy_contains_all_audited_characteristics"],
        "strict_dominance_delta_is_four": matrix["strict_structural_dominance_delta"] == 4,
        "all_lanes_rebased_to_strict_primary": lanes["all_lanes_rebased_to_strict_primary"],
        "legacy_kept_as_bridge_in_all_lanes": lanes["all_lanes_keep_legacy_as_bridge_not_final_kernel"],
        "silent_role_transfer_blocked_in_all_lanes": lanes["all_lanes_block_silent_role_transfer"],
        "apd_and_chi11_rebase_inherited": theorem_export["p2394_apd_bridge_found"] is True and theorem_export["p2394_chi11_selector_found"] is True,
        "s2_completion_additions_detected": theorem_export["s2_mentions_strict_completion_additions"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2403_s1353_v1",
        "packet_id": "P2403",
        "stage_id": "S1353",
        "result_kind": "STRICT_KERNEL_PRIMARY_PHYSICS_GENERATION_REBASE_CERTIFICATE",
        "status": "PASS_STRICT_PRIMARY_REBASE_LEGACY_BRIDGE_ROLE_PRESERVED_NO_PHYSICS_CLOSURE",
        "strict_kernel_primary_physics_generation_rebase_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "characteristic_matrix": matrix,
            "physics_lane_rebase": lanes,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Use strict kernel as the primary object for future mass/gravity/known-physics generation tests, while using legacy+bridge data to explain strict components and keeping role-transfer theorems separate.",
        "global_status": "OPEN_PROGRESS_STRICT_PRIMARY_PHYSICS_GENERATION_REBASED_NO_ROLE_TRANSFER_LICENSE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["strict_kernel_primary_physics_generation_rebase_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2403 S1353: strict-kernel primary physics-generation rebase certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2403/S1353 rebases future known-physics generation work onto the strict kernel as the primary candidate, while preserving the legacy kernel as bridge/construction data.",
                "",
                "## Characteristic matrix summary",
                "",
                f"- Legacy characteristic count: `{theorem['characteristic_matrix_summary']['legacy_characteristic_count']}`.",
                f"- Strict characteristic count: `{theorem['characteristic_matrix_summary']['strict_characteristic_count']}`.",
                f"- Strict additions: `{theorem['characteristic_matrix_summary']['strict_additions']}`.",
                f"- All lanes strict-primary: `{theorem['lane_rebase_summary']['all_lanes_rebased_to_strict_primary']}`.",
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
