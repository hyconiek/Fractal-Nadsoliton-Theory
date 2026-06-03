#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2492_s1442_legacy_strict_claim_specific_minimal_completion_package_certificate.json"
MD = GEN / "p2492_s1442_legacy_strict_claim_specific_minimal_completion_package_certificate.md"

SOURCE_FILES = {
    "P2490_TWO_STAGE_LATTICE": GEN / "p2490_s1440_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate.json",
    "P2491_PIVOTALITY_MATRIX": GEN / "p2491_s1441_legacy_strict_bridge_role_claim_pivotality_matrix_certificate.json",
    "P2434_ROLE_TRANSFER_CLAIM_LATTICE": GEN / "p2434_s1384_conditional_legacy_role_transfer_claim_lattice_certificate.json",
}


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2492|S1442|claim-specific minimal completion package|minimal bridge-role completion package|shared-core completion package|completion package poset|claim package overlap",
        "precursor_packets": "P2490|S1440|two-stage bridge-role closure lattice|P2491|S1441|bridge-role pivotality matrix|P2434|S1384|conditional legacy role-transfer claim lattice",
        "nearby_frontier_cut_language": "theorem-frontier cut|minimal proof package|frontier cut|target-signature lattice|all-open-atom frontier cut",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def stage_counts(package: set[str]) -> dict[str, int]:
    counts = {"bridge": 0, "role": 0}
    for atom in package:
        stage = atom.split("::", 1)[0]
        counts[stage] = counts.get(stage, 0) + 1
    return counts


def classify_role_atom(atom: str, shared_role_atoms: set[str]) -> str:
    role_name = atom.split("::", 1)[1]
    if atom in shared_role_atoms:
        return "post_bridge_common_role_gate"
    if role_name == "alpha_geo_strict_role_successor_theorem":
        return "claim_specific_alpha_geo_successor_gate"
    if role_name == "beta_tors_strict_role_successor_theorem":
        return "claim_specific_beta_tors_successor_gate"
    if role_name == "strict_nonlinear_hierarchy_successor_theorem":
        return "claim_specific_nonlinear_hierarchy_successor_gate"
    if role_name == "chi11_orientation_role_successor_theorem":
        return "claim_specific_chi11_orientation_successor_gate"
    return "claim_specific_role_gate"


def package_certificate(p2490: dict[str, Any], p2491: dict[str, Any], p2434: dict[str, Any]) -> dict[str, Any]:
    lattice = p2490["two_stage_closure_lattice"]
    bridge_atoms = [f"bridge::{atom}" for atom in lattice["bridge_atoms"]]
    role_requirements = p2491["bridge_role_claim_pivotality_matrix"]["claim_requirements"]
    physical_claims = sorted(claim for claim in role_requirements if claim != "all_four_audited_legacy_role_claims")
    packages = {claim: set(role_requirements[claim]) for claim in physical_claims}
    all_roles_package = set(role_requirements["all_four_audited_legacy_role_claims"])
    shared_core = set.intersection(*(packages[claim] for claim in physical_claims))
    union_package = set.union(*(packages[claim] for claim in physical_claims))
    bridge_shared_core = sorted(atom for atom in shared_core if atom.startswith("bridge::"))
    role_shared_core = sorted(atom for atom in shared_core if atom.startswith("role::"))
    physical_package_rows = []
    for claim in physical_claims:
        package = packages[claim]
        role_atoms = sorted(atom for atom in package if atom.startswith("role::"))
        claim_specific = sorted(package - shared_core)
        row = {
            "claim": claim,
            "minimal_current_completion_size": len(package),
            "minimal_current_completion_package": sorted(package),
            "stage_counts": stage_counts(package),
            "shared_core_atom_count": len(shared_core),
            "shared_core_atoms": sorted(shared_core),
            "claim_specific_atom_count": len(claim_specific),
            "claim_specific_atoms": claim_specific,
            "role_atom_classification": [
                {"atom": atom, "class": classify_role_atom(atom, set(role_shared_core))}
                for atom in role_atoms
            ],
        }
        physical_package_rows.append(row)
    physical_package_rows.sort(key=lambda row: (row["minimal_current_completion_size"], row["claim"]))

    overlap_rows = []
    for left, right in combinations(physical_claims, 2):
        left_pkg = packages[left]
        right_pkg = packages[right]
        intersection = left_pkg & right_pkg
        union = left_pkg | right_pkg
        overlap_rows.append({
            "claim_pair": [left, right],
            "intersection_size": len(intersection),
            "union_size": len(union),
            "left_package_is_subset_of_right": left_pkg < right_pkg,
            "right_package_is_subset_of_left": right_pkg < left_pkg,
            "symmetric_difference_atoms": sorted(left_pkg ^ right_pkg),
            "jaccard_fraction": f"{len(intersection)}/{len(union)}",
        })

    package_sizes = [row["minimal_current_completion_size"] for row in physical_package_rows]
    weak_pkg = packages["legacy_weak_mixing_angle"]
    inverse_pkg = packages["legacy_inverse_alpha_em"]
    all_roles_equals_union = all_roles_package == union_package
    role_transfer_ltotal_atoms = {"role::role_transfer_audit_license", "role::role_bearing_ltotal_export"}
    beta_dependent_claims = sorted(
        claim for claim, package in packages.items()
        if "role::beta_tors_strict_role_successor_theorem" in package
    )
    alpha_dependent_claims = sorted(
        claim for claim, package in packages.items()
        if "role::alpha_geo_strict_role_successor_theorem" in package
    )
    return {
        "bridge_atom_count": len(bridge_atoms),
        "physical_claim_count": len(physical_claims),
        "physical_claim_ids": physical_claims,
        "shared_core_atom_count": len(shared_core),
        "shared_core_atoms": sorted(shared_core),
        "bridge_shared_core_atom_count": len(bridge_shared_core),
        "bridge_shared_core_atoms": bridge_shared_core,
        "role_shared_core_atom_count": len(role_shared_core),
        "role_shared_core_atoms": role_shared_core,
        "physical_claim_minimal_completion_rows": physical_package_rows,
        "minimal_physical_claim_completion_size": min(package_sizes),
        "maximal_physical_claim_completion_size": max(package_sizes),
        "all_roles_minimal_completion_size": len(all_roles_package),
        "all_roles_minimal_completion_package": sorted(all_roles_package),
        "union_of_physical_claim_packages_size": len(union_package),
        "union_of_physical_claim_packages": sorted(union_package),
        "all_roles_package_equals_union_of_physical_claim_packages": all_roles_equals_union,
        "package_overlap_rows": overlap_rows,
        "weak_mixing_package_is_proper_subset_of_inverse_alpha_package": weak_pkg < inverse_pkg,
        "role_transfer_and_ltotal_in_every_physical_package": all(role_transfer_ltotal_atoms <= packages[claim] for claim in physical_claims),
        "every_physical_package_contains_all_bridge_atoms": all(set(bridge_atoms) <= packages[claim] for claim in physical_claims),
        "beta_tors_successor_claims": beta_dependent_claims,
        "alpha_geo_successor_claims": alpha_dependent_claims,
        "package_fingerprint_sha256": sha256_json(physical_package_rows),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2492/S1442 legacy-strict claim-specific minimal completion package certificate

`P2492/S1442` converts the P2490 two-stage bridge/role lattice and P2491 claim-pivotality matrix into exact current-state minimal completion packages for each audited legacy physical-role claim.  From the current empty export state, every audited physical claim still requires all eight bridge atoms plus the two common post-bridge role gates `role_transfer_audit_license` and `role_bearing_ltotal_export`.  The weakest package is the weak-mixing successor package of size eleven; the inverse-alpha package strictly extends it by the `beta_tors` successor, while the gravity and torsion/chi11 packages require their own successor atoms.  The all-role package is exactly the union of the four physical-claim packages and has size fourteen.

This is a package/minimality certificate only.  It exports no bridge atom, no selector/source theorem, no QW-2191 discharge, no role-transfer license, no role-bearing `L_total`, no physical-value generator, no legacy-role transfer, and no ToE closure.
"""
    lag_section = """
## P2492/S1442 claim-specific minimal completion package guard

`P2492/S1442` adds a claim-specific minimal-package guard behind `L_total`: any audited legacy-role successor must first carry the full eight-atom bridge package and the shared post-bridge role-transfer/role-bearing-`L_total` core.  Claim-specific atoms then separate weak-mixing, inverse-alpha, gravity-hierarchy, and torsion-to-chi11 successor routes.  The guard is finite bookkeeping/minimality evidence and does not export selector/source closure, QW-2191 discharge, role-transfer, physical-value generation, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2492/S1442 legacy-strict claim-specific minimal completion package certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2492/S1442 claim-specific minimal completion package guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2490 = theorem(sources["P2490_TWO_STAGE_LATTICE"], "legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate")
    p2491 = theorem(sources["P2491_PIVOTALITY_MATRIX"], "legacy_strict_bridge_role_claim_pivotality_matrix_certificate")
    p2434 = theorem(sources["P2434_ROLE_TRANSFER_CLAIM_LATTICE"], "conditional_legacy_role_transfer_claim_lattice_certificate")
    packages = package_certificate(p2490, p2491, p2434)
    theorem_export = {
        "theorem_name": "P2492_T1_legacy_strict_claim_specific_minimal_completion_package_certificate",
        "audited_chain": ["P2490/S1440", "P2491/S1441", "P2434/S1384"],
        "claim_specific_minimal_completion_package_certificate": packages,
        "physical_claim_count": packages["physical_claim_count"],
        "shared_core_atom_count": packages["shared_core_atom_count"],
        "minimal_physical_claim_completion_size": packages["minimal_physical_claim_completion_size"],
        "maximal_physical_claim_completion_size": packages["maximal_physical_claim_completion_size"],
        "all_roles_minimal_completion_size": packages["all_roles_minimal_completion_size"],
        "all_roles_package_equals_union_of_physical_claim_packages": packages["all_roles_package_equals_union_of_physical_claim_packages"],
        "weak_mixing_package_is_proper_subset_of_inverse_alpha_package": packages["weak_mixing_package_is_proper_subset_of_inverse_alpha_package"],
        "role_transfer_and_ltotal_in_every_physical_package": packages["role_transfer_and_ltotal_in_every_physical_package"],
        "every_physical_package_contains_all_bridge_atoms": packages["every_physical_package_contains_all_bridge_atoms"],
        "source_obligation_discharge_exported_by_this_certificate": False,
        "selector_source_theorem_exported_by_this_certificate": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "strict_physical_value_generator_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2492 computes exact current-state minimal completion packages; it does not prove any missing package atom.",
            "The shared ten-atom core is necessary for audited role successors, not sufficient for any physical-role transfer by itself.",
            "The weak-mixing package being smaller than the inverse-alpha package does not license weak-mixing transfer without the bridge and role-transfer theorem.",
            "No QW-2191 discharge, selector/source theorem, role-transfer theorem, physical-value generator, legacy-role transfer, or ToE closure is exported.",
        ],
        "next_honest_step": "Pick an actual theorem-producing target inside the shared core rather than another package-only refinement: either a bridge atom with source content, especially chi11_selector_source_theorem/QW-2191 source if selector work is attempted, or the common post-bridge role_transfer_audit_license/role_bearing_ltotal_export pair.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "physical_claim_count_exact": theorem_export["physical_claim_count"] == 4,
        "shared_core_size_exact": theorem_export["shared_core_atom_count"] == 10,
        "weakest_package_size_exact": theorem_export["minimal_physical_claim_completion_size"] == 11,
        "largest_physical_package_size_exact": theorem_export["maximal_physical_claim_completion_size"] == 12,
        "all_roles_package_size_exact": theorem_export["all_roles_minimal_completion_size"] == 14,
        "all_roles_union_checked": theorem_export["all_roles_package_equals_union_of_physical_claim_packages"],
        "common_gates_checked": theorem_export["role_transfer_and_ltotal_in_every_physical_package"] and theorem_export["every_physical_package_contains_all_bridge_atoms"],
        "known_subset_relation_checked": theorem_export["weak_mixing_package_is_proper_subset_of_inverse_alpha_package"],
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "source_obligation_discharge_exported_by_this_certificate",
            "selector_source_theorem_exported_by_this_certificate",
            "qw2191_discharged_by_this_certificate",
            "role_transfer_licensed_by_this_certificate",
            "role_bearing_ltotal_exported",
            "strict_physical_value_generator_exported",
            "legacy_role_transfer_exported",
            "toe_closure_exported",
        ]),
    }
    return {
        "packet_id": "P2492",
        "stage_id": "S1442",
        "status": "LEGACY_STRICT_CLAIM_SPECIFIC_MINIMAL_COMPLETION_PACKAGE_CERTIFICATE_NO_ATOM_EXPORT_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_strict_claim_specific_minimal_completion_package_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_strict_claim_specific_minimal_completion_package_certificate"]["theorem_export"]
    cert = t["claim_specific_minimal_completion_package_certificate"]
    lines = [
        "# P2492/S1442 legacy-strict claim-specific minimal completion package certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Minimal package summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Physical claims: `{t['physical_claim_count']}`.",
        f"Shared current-state core atoms: `{t['shared_core_atom_count']}`.",
        f"Minimal physical-claim completion size: `{t['minimal_physical_claim_completion_size']}`.",
        f"Maximal physical-claim completion size: `{t['maximal_physical_claim_completion_size']}`.",
        f"All-role completion size: `{t['all_roles_minimal_completion_size']}`.",
        f"All-role package equals union of physical packages: `{t['all_roles_package_equals_union_of_physical_claim_packages']}`.",
        "",
        "## Claim rows",
        "",
    ]
    for row in cert["physical_claim_minimal_completion_rows"]:
        lines.append(
            f"- `{row['claim']}`: size `{row['minimal_current_completion_size']}`, "
            f"stage counts `{row['stage_counts']}`, claim-specific atoms `{', '.join(row['claim_specific_atoms'])}`."
        )
    lines += [
        "",
        "## Shared core",
        "",
        f"Bridge shared core count: `{cert['bridge_shared_core_atom_count']}`.",
        f"Role shared core atoms: `{', '.join(cert['role_shared_core_atoms'])}`.",
        "",
        "## Negative controls",
        "",
        "P2492 does not export a bridge atom, selector/source theorem, QW-2191 discharge, role-transfer license, role-bearing L_total, physical-value generator, legacy-role transfer, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2492 computes the smallest exact packages that would have to be supplied before each audited legacy role claim could even be considered transferable.  It narrows proof planning without pretending that any missing proof has been supplied.",
        "",
        "## Fingerprints",
        "",
        f"Package fingerprint: `{cert['package_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['legacy_strict_claim_specific_minimal_completion_package_certificate']['theorem_fingerprint_sha256']}`.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
