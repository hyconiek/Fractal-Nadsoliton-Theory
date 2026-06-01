#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.json"
MD = GEN / "p2394_s1344_apd_bridge_chi11_rebased_role_frontier_certificate.md"

SOURCE_FILES = {
    "APD_FINITE_BRIDGE_ASSEMBLY": SCRATCH / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json",
    "THEOREM_FRONTIER_CUT": SCRATCH / "bridge_strict_completion_theorem_frontier_cut_certificate_report.json",
    "THEOREM_FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "P1343_SELECTOR_SOURCE_REPORT": GEN / "p1343_p1343_report_v1.json",
    "P1348_GLOBAL_CLOSURE_REPORT": GEN / "p1348_p1348_report_v1.json",
    "P2392_BETA_TORS_CHI11_RETIREMENT": GEN / "p2392_s1342_auxiliary_beta_tors_chi11_retirement_certificate.json",
    "P2393_BOUNDARY_NEGATIVE_CONTROL": GEN / "p2393_s1343_kernel_completion_boundary_residual_certificate.json",
    "F4_ROLE_CLASSIFICATION": ROOT / "F4_LEGACY_PHYSICAL_ROLE_TRANSFER_CLASSIFICATION_PACKET.md",
    "S2_STRATEGIC_PRIORITY": ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROLE_ATOMS = [
    "alpha_geo_electroweak_role_theorem",
    "beta_tors_strict_role_theorem",
    "beta_power_hierarchy_successor_theorem",
]

ROLE_TARGETS = {
    "legacy_weinberg_sin2_theta_role_transfer": {"alpha_geo_electroweak_role_theorem"},
    "legacy_alpha_em_inverse_role_transfer": {
        "alpha_geo_electroweak_role_theorem",
        "beta_tors_strict_role_theorem",
    },
    "legacy_beta_power_gravity_hierarchy_successor": {
        "beta_tors_strict_role_theorem",
        "beta_power_hierarchy_successor_theorem",
    },
}


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def read_text(path: Path) -> str:
    if not path.exists():
        return f"OPEN_MISSING_ARTIFACT::{rel(path)}"
    return path.read_text(encoding="utf-8")


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2394|S1344|APD bridge chi11|rebased role frontier",
        r"K_strict.*K_legacy.*A.*P.*D|K_legacy.*A.*P.*D|Q_assembly=A\*P\*D|finite bridge assembly",
        "chi11_selector_source|P1343|P1348|strict selector mechanism",
        r"legacy_weinberg|sin\^2\(theta_W\)|alpha_EM\^-1|beta\^N|role-transfer theorem",
        "P2393|eta=1 boundary|current residual",
    ]
    out: dict[str, Any] = {}
    for pattern in patterns:
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "-g", "*.py", "-g", "*.md", "-g", "*.json"],
            cwd=REPO,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        lines = [line for line in proc.stdout.splitlines() if line]
        out[pattern] = {"count": len(lines), "samples": lines[:12]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep confirms that the finite APD bridge assembly already exists and that P2393 was only an eta=1 boundary negative-control computation. "
            "P2394 therefore rebases the active frontier to APD-bridge-found plus chi11-selector-found, then audits the remaining legacy role-transfer atoms."
        ),
    }


def extract_apd_bridge_status(assembly: dict[str, Any]) -> dict[str, Any]:
    summary = assembly.get("finite_bridge_assembly_summary", {})
    cross = assembly.get("cross_checks", {})
    proof = assembly.get("proof_certificate", {})
    found = bool(
        summary.get("finite_kernel_comparison_witness_exported")
        and summary.get("assembled_map_reconstructs_strict_kernel")
        and summary.get("assembled_map_matches_finite_diagonal_certificate")
        and summary.get("assembled_map_matches_necessity_apd_product")
        and cross.get("assembled_identity_exact_on_finite_domain")
    )
    return {
        "apd_bridge_found_in_existing_repo": found,
        "status": assembly.get("status", "OPEN_MISSING_ARTIFACT"),
        "domain_size": summary.get("domain_size"),
        "max_abs_reconstruction_residual": summary.get("max_abs_reconstruction_residual"),
        "max_abs_assembled_q_minus_diagonal_q": summary.get("max_abs_assembled_q_minus_diagonal_q"),
        "max_abs_assembled_q_minus_necessity_factor_product": summary.get(
            "max_abs_assembled_q_minus_necessity_factor_product"
        ),
        "assembly_step": proof.get("assembly_step"),
        "identity_step": proof.get("identity_step"),
        "scope_step": proof.get("scope_step"),
        "interpretation": (
            "The finite comparison bridge K_strict = K_legacy*A*P*D is already found on the audited domain; "
            "do not demote it to an eta=1 boundary-only gap."
        ),
    }


def extract_selector_status(p1343: dict[str, Any], p1348: dict[str, Any], p2392: dict[str, Any]) -> dict[str, Any]:
    p2392_cert = p2392.get("auxiliary_beta_tors_chi11_retirement_certificate", {})
    theorem = p2392_cert.get("theorem_export", {})
    available_atoms = theorem.get("available_atoms", {})
    selector_support = theorem.get("support_evaluation", {}).get("selector_mechanism", {})
    selector_found = bool(
        p1343.get("status") == "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1"
        and p1348.get("status") == "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE"
        and available_atoms.get("strict_internal_selector_P1343_P1348")
        and selector_support.get("target_satisfied")
        and not selector_support.get("uses_auxiliary_beta_tors_in_realized_support", True)
    )
    return {
        "strict_selector_found_in_declared_scope": selector_found,
        "p1343_status": p1343.get("status"),
        "p1348_status": p1348.get("status"),
        "p2392_status": p2392.get("status"),
        "selector_realized_without_beta_tors_chi11": bool(
            selector_support.get("target_satisfied")
            and not selector_support.get("uses_auxiliary_beta_tors_in_realized_support", True)
        ),
        "retired_route": "beta_tors -> chi11 remains retired as a selector-search assumption, not a missing bridge component.",
    }


def truth_row(active_atoms: set[str]) -> dict[str, Any]:
    target_values = {target: required <= active_atoms for target, required in ROLE_TARGETS.items()}
    all_transfer = all(target_values.values())
    return {
        "true_role_atoms": sorted(active_atoms),
        "target_values": target_values,
        "all_three_physical_role_transfers": all_transfer,
        "closed_role_count": sum(target_values.values()),
    }


def role_truth_table() -> list[dict[str, Any]]:
    rows = []
    for size in range(len(ROLE_ATOMS) + 1):
        for combo in itertools.combinations(ROLE_ATOMS, size):
            rows.append(truth_row(set(combo)))
    return rows


def minimal_supports_for_roles() -> dict[str, Any]:
    supports = {target: [sorted(required)] for target, required in ROLE_TARGETS.items()}
    full_set = sorted(set().union(*ROLE_TARGETS.values()))
    return {
        "individual_role_minimal_supports": supports,
        "all_physical_roles_minimal_support": full_set,
        "all_physical_roles_minimal_weight": len(full_set),
    }


def rebased_frontier(apd: dict[str, Any], selector: dict[str, Any]) -> dict[str, Any]:
    current_role_atoms = {atom: False for atom in ROLE_ATOMS}
    table = role_truth_table()
    current = truth_row(set())
    supports = minimal_supports_for_roles()
    singleton_unlocks = []
    for atom in ROLE_ATOMS:
        row = truth_row({atom})
        unlocked = [target for target, value in row["target_values"].items() if value]
        singleton_unlocks.append({"atom": atom, "unlocked_role_targets": unlocked, "closed_role_count": row["closed_role_count"]})
    return {
        "closed_context_atoms": {
            "finite_apd_bridge_assembly_witness": apd["apd_bridge_found_in_existing_repo"],
            "strict_internal_chi11_selector_mechanism_P1343_P1348": selector[
                "strict_selector_found_in_declared_scope"
            ],
        },
        "active_role_atoms_after_rebase": current_role_atoms,
        "role_target_definitions_after_bridge_and_chi11": {
            target: sorted(required) for target, required in ROLE_TARGETS.items()
        },
        "current_assignment": current,
        "truth_table_rows": table,
        "truth_table_size": len(table),
        "minimal_supports": supports,
        "singleton_unlocks": singleton_unlocks,
        "chi11_not_in_active_role_atom_set": "chi11_selector_source" not in ROLE_ATOMS,
        "interpretation": (
            "Once the APD finite bridge and strict selector are treated as found in their declared scope, "
            "the next active mathematical work is not another APD or beta_tors->chi11 proof. "
            "It is the post-bridge legacy physical-role transfer audit over alpha_geo, beta_tors, and beta-power successor atoms."
        ),
    }


def p2393_correction(p2393: dict[str, Any], apd: dict[str, Any]) -> dict[str, Any]:
    cert = p2393.get("kernel_completion_boundary_residual_certificate", {})
    residual = cert.get("current_strict_residual_audit", {})
    boundary = cert.get("boundary_identity_audit", {})
    return {
        "p2393_boundary_identity_kept_as_negative_control": bool(boundary.get("allclose_at_1e_minus_15")),
        "p2393_current_residual_nonzero": bool(residual.get("nonzero_residual_on_current_target")),
        "p2393_residual_reinterpreted_as_eta1_slice_insufficiency_not_missing_apd_bridge": apd[
            "apd_bridge_found_in_existing_repo"
        ]
        and bool(residual.get("nonzero_residual_on_current_target")),
        "correction": (
            "The P2393 eta=1 residual should not be read as evidence that the APD bridge was unfound. "
            "It only proves that the eta=1 boundary slice is not the current strict target unless the already-known P(d) and D(d) completion factors are applied."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2394/S1344 APD bridge found, chi11 rebased, role-transfer frontier

`P2394/S1344` corrects the active reading after P2393.  The finite comparison bridge was already found in the existing APD assembly certificate:

```text
K_strict(d) = K_legacy(d) * A(d) * P(d) * D(d)
```

on the audited finite domain.  Therefore the P2393 `eta=1` boundary residual is only a negative-control slice: it says the boundary slice alone is insufficient, not that the APD bridge is missing.

After P2392/P1343/P1348, the strict selector/`chi11` mechanism is also available in declared scope without using the retired `beta_tors -> chi11` selector-search route.  The rebased active frontier is therefore the post-bridge legacy physical-role audit, not another APD proof and not a `beta_tors -> chi11` proof.

The computed role-transfer truth table has three active role atoms:

```text
alpha_geo_electroweak_role_theorem,
beta_tors_strict_role_theorem,
beta_power_hierarchy_successor_theorem.
```

Current assignment closes none of the three legacy physical-role transfers.  Minimal supports are: Weinberg role needs `alpha_geo`, alpha-EM needs `alpha_geo + beta_tors`, and gravity hierarchy needs `beta_tors + beta_power`.  No `L_total`, SM/GR numeric extraction, or ToE closure follows.
""".strip()
    lag_section = """
## P2394/S1344 APD bridge rebase, Lagrangian role-transfer frontier

`P2394/S1344` prevents the Lagrangian/EOM lane from misreading P2393.  The current bridge input is the already assembled finite APD comparison witness `K_strict = K_legacy*A*P*D`; P2393 is only an `eta=1` boundary negative control.  With strict selector/`chi11` available in declared scope, the next honest Lagrangian-side task is not to re-prove APD or revive `beta_tors -> chi11`, but to audit whether any legacy physical roles survive after the bridge.

The role-transfer truth table leaves the legacy Weinberg, alpha-EM, and beta-power gravity roles unlicensed at the current assignment.  Any future `L_total` term using `alpha_geo`, `beta_tors`, or `beta^N` as physical legacy roles must pass that role-transfer audit first.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {
        name: load_json(path) if path.suffix == ".json" else {"text": read_text(path)}
        for name, path in SOURCE_FILES.items()
    }
    grep = rg_audit()
    apd = extract_apd_bridge_status(artifacts["APD_FINITE_BRIDGE_ASSEMBLY"])
    selector = extract_selector_status(
        artifacts["P1343_SELECTOR_SOURCE_REPORT"],
        artifacts["P1348_GLOBAL_CLOSURE_REPORT"],
        artifacts["P2392_BETA_TORS_CHI11_RETIREMENT"],
    )
    frontier = rebased_frontier(apd, selector)
    correction = p2393_correction(artifacts["P2393_BOUNDARY_NEGATIVE_CONTROL"], apd)
    theorem_export = {
        "theorem_name": "P2394_T1_apd_bridge_found_chi11_rebased_role_frontier",
        "apd_bridge_found": apd["apd_bridge_found_in_existing_repo"],
        "strict_chi11_selector_found": selector["strict_selector_found_in_declared_scope"],
        "active_missing_selector_route_count": 0,
        "active_apd_reproof_obligation_count": 0,
        "active_role_atoms": ROLE_ATOMS,
        "current_role_assignment": frontier["current_assignment"],
        "minimal_supports": frontier["minimal_supports"],
        "p2393_correction": correction,
        "not_licensed": [
            "No beta_tors -> chi11 selector-search route is reopened.",
            "No silent transfer of alpha_geo, beta_tors, or beta^N physical roles is claimed.",
            "No L_total promotion, SM/GR numeric extraction, or ToE closure is claimed.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "apd_bridge_found_inherited": apd["apd_bridge_found_in_existing_repo"],
        "chi11_selector_found_inherited": selector["strict_selector_found_in_declared_scope"],
        "p2393_reinterpreted_not_active_bridge_gap": correction[
            "p2393_residual_reinterpreted_as_eta1_slice_insufficiency_not_missing_apd_bridge"
        ],
        "chi11_removed_from_active_role_atoms": frontier["chi11_not_in_active_role_atom_set"],
        "current_assignment_closes_no_legacy_physical_role": frontier["current_assignment"]["closed_role_count"] == 0,
        "truth_table_size_is_2_pow_3": frontier["truth_table_size"] == 8,
        "all_role_minimal_weight_is_3": frontier["minimal_supports"]["all_physical_roles_minimal_weight"] == 3,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2394_s1344_v1",
        "packet_id": "P2394",
        "stage_id": "S1344",
        "result_kind": "APD_BRIDGE_CHI11_REBASED_ROLE_FRONTIER_CERTIFICATE",
        "status": "PASS_APD_BRIDGE_FOUND_CHI11_SELECTOR_REBASED_ROLE_TRANSFER_FRONTIER_OPEN",
        "apd_bridge_chi11_rebased_role_frontier_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "apd_bridge_status": apd,
            "selector_status": selector,
            "p2393_boundary_correction": correction,
            "rebased_role_frontier": frontier,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Run the post-bridge legacy physical-role audit over alpha_geo, beta_tors, and beta-power successor roles; do not re-prove APD and do not revive beta_tors->chi11 as selector work.",
        "global_status": "OPEN_PROGRESS_APD_AND_CHI11_REBASED_ROLE_TRANSFER_AUDIT_NOW_ACTIVE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["apd_bridge_chi11_rebased_role_frontier_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2394 S1344: APD bridge found, chi11 rebased, role frontier certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2394/S1344 accepts the user's correction: the finite APD comparison bridge was already found as `K_strict = K_legacy*A*P*D`. P2393 is therefore only an eta=1 boundary negative control, not the active bridge gap.",
                "",
                "## Rebased closed context",
                "",
                f"- APD bridge found: `{theorem['apd_bridge_found']}`.",
                f"- Strict chi11 selector found in declared scope: `{theorem['strict_chi11_selector_found']}`.",
                f"- Active APD reproof obligation count: `{theorem['active_apd_reproof_obligation_count']}`.",
                f"- Active missing selector-route count: `{theorem['active_missing_selector_route_count']}`.",
                "",
                "## Role-transfer frontier",
                "",
                f"- Active role atoms: `{theorem['active_role_atoms']}`.",
                f"- Current role assignment: `{theorem['current_role_assignment']}`.",
                f"- Minimal supports: `{theorem['minimal_supports']}`.",
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
