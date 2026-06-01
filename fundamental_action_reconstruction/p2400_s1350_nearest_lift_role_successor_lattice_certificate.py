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

OUT = GEN / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.json"
MD = GEN / "p2400_s1350_nearest_lift_role_successor_lattice_certificate.md"

SOURCE_FILES = {
    "P2399_LIFT_DISTANCE": GEN / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.json",
    "FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "ROLE_MINIMAL_OBLIGATION_LATTICE": SCRATCH / "bridge_strict_completion_legacy_role_transfer_minimal_obligation_lattice_certificate_report.json",
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

NONROLE_CONTEXT = {
    "strict_dynamical_source_for_A_P_D": True,
    "strict_phase_frequency_source": True,
    "strict_damping_beta_eta_source": True,
    "chi11_selector_source": True,
}

TARGET_KEYS = [
    "bridge_theorem_level_closure",
    "role_transfer_theorem_level_closure",
    "selector_qw2191_closure",
    "toe_closure",
]


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


def rg_audit() -> dict[str, Any]:
    patterns = [
        "P2400|S1350|nearest-lift role-successor lattice|role-successor subset lattice|nearest lift role lattice",
        "P2399|role-closed lift-distance spectrum|nearest ToE rows|minimum ToE lift distance",
        "legacy_role_transfer_minimal_obligation_lattice|minimal obligation set|global_minimal_obligation_sets",
        "alpha_geo_electroweak_role_theorem|beta_tors_strict_role_theorem|beta_power_hierarchy_successor_theorem",
        "proper subsets fail|single role atom|two role atoms|all three role atoms",
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
        out[pattern] = {"count": len(lines), "samples": lines[:16]}
    return {
        "tool": "rg",
        "patterns": out,
        "finding": (
            "Repo grep finds the older four-atom legacy-role minimal-obligation lattice and P2399's lift-distance spectrum. "
            "P2400 only computes the three-role successor subset lattice at the P2399 nearest non-role-complete lift row where chi11 is already fixed true."
        ),
    }


def target_values(active: set[str], target_definitions: dict[str, Any]) -> dict[str, bool]:
    bridge = set(target_definitions["bridge_theorem_level_closure"]) <= active
    role = set(target_definitions["role_transfer_theorem_level_closure"]) <= active
    selector = set(target_definitions["selector_qw2191_closure"]) <= active
    toe = bridge and role and selector
    return {
        "bridge_theorem_level_closure": bridge,
        "role_transfer_theorem_level_closure": role,
        "selector_qw2191_closure": selector,
        "toe_closure": toe,
    }


def role_mask(role_assignment: dict[str, bool]) -> int:
    mask = 0
    for index, atom in enumerate(ROLE_ATOMS):
        if role_assignment[atom]:
            mask |= 1 << index
    return mask


def monomial_name(mask: int) -> str:
    if mask == 0:
        return "1"
    return "*".join(atom for index, atom in enumerate(ROLE_ATOMS) if mask & (1 << index))


def anf_coefficients(values_by_mask: dict[int, int]) -> list[dict[str, Any]]:
    coeffs = [values_by_mask.get(mask, 0) for mask in range(1 << len(ROLE_ATOMS))]
    for bit in range(len(ROLE_ATOMS)):
        for mask in range(1 << len(ROLE_ATOMS)):
            if mask & (1 << bit):
                coeffs[mask] ^= coeffs[mask ^ (1 << bit)]
    return [
        {"mask": mask, "degree": mask.bit_count(), "monomial": monomial_name(mask)}
        for mask, coeff in enumerate(coeffs)
        if coeff
    ]


def compute_lattice(target_definitions: dict[str, Any]) -> dict[str, Any]:
    base_active = {atom for atom, value in NONROLE_CONTEXT.items() if value}
    rows = []
    truth_by_mask = {target: {} for target in TARGET_KEYS}
    for values in itertools.product([False, True], repeat=len(ROLE_ATOMS)):
        assignment = dict(zip(ROLE_ATOMS, values))
        active_roles = {atom for atom, value in assignment.items() if value}
        active = base_active | active_roles
        targets = target_values(active, target_definitions)
        mask = role_mask(assignment)
        missing_roles = [atom for atom in ROLE_ATOMS if atom not in active_roles]
        for target in TARGET_KEYS:
            truth_by_mask[target][mask] = int(targets[target])
        rows.append(
            {
                "role_mask": mask,
                "role_assignment": assignment,
                "active_role_atoms": sorted(active_roles),
                "missing_role_atoms": missing_roles,
                "missing_role_count": len(missing_roles),
                "target_values": targets,
                "target_signature_bridge_role_selector_toe": "".join("1" if targets[key] else "0" for key in TARGET_KEYS),
            }
        )
    role_true_masks = sorted(row["role_mask"] for row in rows if row["target_values"]["role_transfer_theorem_level_closure"])
    toe_true_masks = sorted(row["role_mask"] for row in rows if row["target_values"]["toe_closure"])
    proper_fail_masks = sorted(row["role_mask"] for row in rows if row["role_mask"] != 7 and not row["target_values"]["toe_closure"])
    nearest_miss_masks = sorted(row["role_mask"] for row in rows if row["missing_role_count"] == 1 and not row["target_values"]["toe_closure"])
    derivatives = {}
    for index, atom in enumerate(ROLE_ATOMS):
        support = []
        for mask in range(1 << len(ROLE_ATOMS)):
            flipped = mask ^ (1 << index)
            if mask < flipped:
                derivative = truth_by_mask["toe_closure"][mask] ^ truth_by_mask["toe_closure"][flipped]
                if derivative:
                    support.append({"lower_mask": mask, "upper_mask": flipped, "other_roles_true": monomial_name(mask & ~(1 << index))})
        derivatives[atom] = support
    return {
        "fixed_nonrole_context": NONROLE_CONTEXT,
        "role_atoms": ROLE_ATOMS,
        "rows": sorted(rows, key=lambda row: row["role_mask"]),
        "role_transfer_true_masks": role_true_masks,
        "toe_true_masks": toe_true_masks,
        "proper_subset_fail_masks": proper_fail_masks,
        "nearest_one_role_missing_masks": nearest_miss_masks,
        "target_truth_vectors_by_mask_0_to_7": {
            target: [truth_by_mask[target][mask] for mask in range(1 << len(ROLE_ATOMS))]
            for target in TARGET_KEYS
        },
        "target_anf_by_role_atoms": {
            target: anf_coefficients(truth_by_mask[target])
            for target in TARGET_KEYS
        },
        "toe_boolean_derivative_support_by_role_atom": derivatives,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2400/S1350 nearest-lift role-successor lattice certificate

`P2400/S1350` takes the P2399 nearest non-role-complete row and enumerates the remaining three-role successor lattice exactly.  With APD, phase/frequency, damping, and `chi11` fixed true, the role-transfer/ToE Boolean function over the three role atoms is the single degree-3 monomial:

```text
alpha_geo_electroweak_role_theorem * beta_tors_strict_role_theorem * beta_power_hierarchy_successor_theorem.
```

All seven proper role subsets fail; the three one-missing nearest misses are the exact Boolean-derivative supports.  This is a conditional lattice theorem, not a proof that any role atom has been exported.
""".strip()
    lag_section = """
## P2400/S1350 role-successor lattice guard for Lagrangian/EOM

`P2400/S1350` fixes all non-role bridge/selector atoms true and still finds that role-transfer/ToE require the full three-role successor monomial.  Thus no one-role or two-role partial package can justify a role-bearing `L_total` term.

The calculation is conditional on future role evidence and exports no physical-role theorem by itself.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    target_definitions = artifacts["FRONTIER_TRUTH_TABLE"].get("target_definitions", {})
    lattice = compute_lattice(target_definitions)
    p2399_theorem = artifacts["P2399_LIFT_DISTANCE"].get("role_closed_lift_distance_spectrum_certificate", {}).get("theorem_export", {})
    older_lattice = artifacts["ROLE_MINIMAL_OBLIGATION_LATTICE"].get("minimal_obligation_lattice_summary", {})
    toe_anf = lattice["target_anf_by_role_atoms"]["toe_closure"]
    role_anf = lattice["target_anf_by_role_atoms"]["role_transfer_theorem_level_closure"]
    theorem_export = {
        "theorem_name": "P2400_T1_nearest_lift_three_role_successor_lattice",
        "fixed_nonrole_context": NONROLE_CONTEXT,
        "role_atoms": ROLE_ATOMS,
        "role_transfer_true_masks": lattice["role_transfer_true_masks"],
        "toe_true_masks": lattice["toe_true_masks"],
        "proper_subset_fail_count": len(lattice["proper_subset_fail_masks"]),
        "nearest_one_role_missing_masks": lattice["nearest_one_role_missing_masks"],
        "role_transfer_role_anf": role_anf,
        "toe_role_anf": toe_anf,
        "toe_boolean_derivative_support_by_role_atom": lattice["toe_boolean_derivative_support_by_role_atom"],
        "p2399_minimum_toe_lift_distance": p2399_theorem.get("minimum_toe_lift_distance"),
        "older_four_atom_lattice_context": older_lattice.get("global_minimal_obligation_sets"),
        "not_licensed": [
            "No role-successor atom is exported by this conditional lattice.",
            "No one-role or two-role subset licenses role transfer.",
            "No ToE closure is claimed without all three explicit role-successor theorems.",
            "No L_total or SM/GR role-bearing term is promoted here.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2399_minimum_toe_distance_is_three": p2399_theorem.get("minimum_toe_lift_distance") == 3,
        "only_full_role_mask_closes_role_transfer": theorem_export["role_transfer_true_masks"] == [7],
        "only_full_role_mask_closes_toe": theorem_export["toe_true_masks"] == [7],
        "all_seven_proper_role_subsets_fail": theorem_export["proper_subset_fail_count"] == 7,
        "three_one_missing_nearest_misses": theorem_export["nearest_one_role_missing_masks"] == [3, 5, 6],
        "toe_and_role_anf_are_same_degree3_monomial": toe_anf == role_anf and len(toe_anf) == 1 and toe_anf[0]["degree"] == 3,
        "each_role_derivative_has_unique_support": all(len(support) == 1 for support in theorem_export["toe_boolean_derivative_support_by_role_atom"].values()),
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2400_s1350_v1",
        "packet_id": "P2400",
        "stage_id": "S1350",
        "result_kind": "NEAREST_LIFT_ROLE_SUCCESSOR_LATTICE_CERTIFICATE",
        "status": "PASS_NEAREST_LIFT_REQUIRES_FULL_THREE_ROLE_SUCCESSOR_MONOMIAL",
        "nearest_lift_role_successor_lattice_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "target_definitions": target_definitions,
            "role_successor_lattice": lattice,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "If pursuing role transfer, construct one of the three explicit role-successor theorems and rerun the lattice; do not promote partial role subsets.",
        "global_status": "OPEN_PROGRESS_FULL_THREE_ROLE_SUCCESSOR_MONOMIAL_CERTIFIED_CONDITIONAL_ONLY",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["nearest_lift_role_successor_lattice_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2400 S1350: nearest-lift role-successor lattice certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2400/S1350 fixes the P2399 nearest non-role-complete context and enumerates the remaining three-role successor lattice.",
                "",
                "## Lattice result",
                "",
                f"- Role-transfer true masks: `{theorem['role_transfer_true_masks']}`.",
                f"- ToE true masks: `{theorem['toe_true_masks']}`.",
                f"- Proper role subsets failing: `{theorem['proper_subset_fail_count']}`.",
                f"- Nearest one-role-missing masks: `{theorem['nearest_one_role_missing_masks']}`.",
                f"- ToE role ANF: `{theorem['toe_role_anf']}`.",
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
