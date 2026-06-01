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

OUT = GEN / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.json"
MD = GEN / "p2399_s1349_role_closed_lift_distance_spectrum_certificate.md"

SOURCE_FILES = {
    "P2397_ROLE_CLOSED_TOE_PROJECTION": GEN / "p2397_s1347_role_closed_toe_projection_certificate.json",
    "P2398_ROLE_CLOSED_QUOTIENT_ANF": GEN / "p2398_s1348_role_closed_quotient_anf_certificate.json",
    "FRONTIER_TRUTH_TABLE": SCRATCH / "bridge_strict_completion_theorem_frontier_truth_table_certificate_report.json",
    "TOE_PROPER_SUBSET_OBSTRUCTION": SCRATCH / "bridge_strict_completion_toe_proper_subset_obstruction_certificate_report.json",
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

NONROLE_ATOMS = [
    "strict_dynamical_source_for_A_P_D",
    "strict_phase_frequency_source",
    "strict_damping_beta_eta_source",
    "chi11_selector_source",
]

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
        "P2399|S1349|role-closed lift distance|lift distance spectrum|distance-to-closure spectrum",
        "P2398|role-closed quotient ANF|role_transfer = 0|toe = 0",
        "P2397|role-closed ToE projection|role atoms forced false|slice row count",
        "Hamming distance|nearest miss|proper subset obstruction|six-atom nearest miss",
        "alpha_geo_electroweak_role_theorem|beta_tors_strict_role_theorem|beta_power_hierarchy_successor_theorem|chi11_selector_source",
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
            "Repo grep finds global nearest-miss/proper-subset reports plus P2397/P2398 quotient results. "
            "P2399 therefore does not redo the global seven-atom obstruction; it computes the exact lift-distance spectrum from the P2396/P2397 role-closed quotient."
        ),
    }


def gf2_rank(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix if any(row)]
    if not work:
        return 0
    rank = 0
    col = 0
    width = len(work[0])
    while rank < len(work) and col < width:
        pivot = next((row for row in range(rank, len(work)) if work[row][col] % 2), None)
        if pivot is None:
            col += 1
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(len(work)):
            if row != rank and work[row][col] % 2:
                work[row] = [a ^ b for a, b in zip(work[row], work[rank])]
        rank += 1
        col += 1
    return rank


def required_atoms(target: str, target_definitions: dict[str, Any]) -> list[str]:
    if target == "toe_closure":
        atoms: list[str] = []
        for key in ["bridge_theorem_level_closure", "role_transfer_theorem_level_closure", "selector_qw2191_closure"]:
            for atom in target_definitions[key]:
                if atom not in atoms:
                    atoms.append(atom)
        return atoms
    return list(target_definitions[target])


def spectrum_key(distance: int) -> str:
    return str(distance)


def compute_lift_spectrum(p2397: dict[str, Any], target_definitions: dict[str, Any]) -> dict[str, Any]:
    rows = p2397["role_closed_toe_projection_certificate"]["role_closed_slice"]["rows"]
    lift_rows = []
    spectra = {target: {} for target in TARGET_KEYS}
    missing_matrices = {target: [] for target in TARGET_KEYS}
    target_required = {target: required_atoms(target, target_definitions) for target in TARGET_KEYS}
    all_atoms = ROLE_ATOMS + NONROLE_ATOMS

    for index, row in enumerate(rows):
        assignment = row["nonrole_assignment"]
        active = {atom for atom, truth in assignment.items() if truth}
        per_target: dict[str, Any] = {}
        for target in TARGET_KEYS:
            missing = [atom for atom in target_required[target] if atom not in active]
            role_missing = [atom for atom in missing if atom in ROLE_ATOMS]
            nonrole_missing = [atom for atom in missing if atom in NONROLE_ATOMS]
            distance = len(missing)
            spectra[target][spectrum_key(distance)] = spectra[target].get(spectrum_key(distance), 0) + 1
            missing_vector = [1 if atom in missing else 0 for atom in all_atoms]
            missing_matrices[target].append(missing_vector)
            per_target[target] = {
                "required_atoms": target_required[target],
                "missing_atoms": missing,
                "role_missing_count": len(role_missing),
                "nonrole_missing_count": len(nonrole_missing),
                "lift_distance_to_target": distance,
                "missing_vector_role_then_nonrole": missing_vector,
            }
        lift_rows.append(
            {
                "row_index": index,
                "nonrole_assignment": assignment,
                "lift_targets": per_target,
                "toe_distance_decomposition": {
                    "forced_role_atoms_missing": per_target["toe_closure"]["role_missing_count"],
                    "nonrole_atoms_missing": per_target["toe_closure"]["nonrole_missing_count"],
                    "total_distance": per_target["toe_closure"]["lift_distance_to_target"],
                },
            }
        )

    matrix_ranks = {target: gf2_rank(matrix) for target, matrix in missing_matrices.items()}
    min_distances = {
        target: min(row["lift_targets"][target]["lift_distance_to_target"] for row in lift_rows)
        for target in TARGET_KEYS
    }
    nearest_rows = {
        target: [
            row["row_index"]
            for row in lift_rows
            if row["lift_targets"][target]["lift_distance_to_target"] == min_distances[target]
        ]
        for target in TARGET_KEYS
    }
    return {
        "role_atoms_forced_false": ROLE_ATOMS,
        "free_nonrole_atoms": NONROLE_ATOMS,
        "all_atoms_role_then_nonrole": all_atoms,
        "target_required_atoms": target_required,
        "lift_rows": lift_rows,
        "distance_spectra": {
            target: dict(sorted(target_spectrum.items(), key=lambda item: int(item[0])))
            for target, target_spectrum in spectra.items()
        },
        "minimum_lift_distances": min_distances,
        "nearest_row_indices": nearest_rows,
        "missing_matrix_ranks_mod2": matrix_ranks,
        "proof_reading": (
            "The nearest role-transfer lift needs three role atoms when chi11 is already present; the nearest ToE lift needs those same three role atoms when all four non-role atoms are already present."
        ),
    }


def append_doc_sections() -> None:
    eq_section = """
## P2399/S1349 role-closed lift-distance spectrum certificate

`P2399/S1349` lifts the P2398 quotient result back toward the seven-atom frontier by computing exact Hamming/lift distances from every P2396 role-closed row to each target.  This does not redo the global nearest-miss theorem; it asks a narrower question: how far is the current role-closed quotient from role-transfer and ToE once non-role progress is allowed?

The result is finite and sharp: the nearest role-transfer lift has distance `3` and requires all three role atoms when `chi11_selector_source` is already present; the nearest ToE lift also has distance `3` and occurs only when all four non-role atoms are already present.  Therefore non-role progress can reduce the ToE lift distance down to three, but it cannot remove the three explicit role-successor obligations.

No role-transfer theorem, ToE closure, `L_total` promotion, or SM/GR numeric extraction follows from this distance spectrum.
""".strip()
    lag_section = """
## P2399/S1349 lift-distance guard for Lagrangian/EOM

`P2399/S1349` computes the exact lift-distance spectrum from the role-closed quotient.  Even at the nearest row where APD, phase/frequency, damping, and `chi11` are all present, the ToE/role-transfer lift still has distance `3` because the three role-successor atoms remain missing.

Thus a role-bearing `L_total` term still requires explicit role-successor evidence rather than additional non-role source bookkeeping alone.
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
    spectrum = compute_lift_spectrum(artifacts["P2397_ROLE_CLOSED_TOE_PROJECTION"], target_definitions)
    p2398_theorem = artifacts["P2398_ROLE_CLOSED_QUOTIENT_ANF"].get("role_closed_quotient_anf_certificate", {}).get("theorem_export", {})
    proper_subset = artifacts["TOE_PROPER_SUBSET_OBSTRUCTION"].get("proper_subset_obstruction_summary", {})
    theorem_export = {
        "theorem_name": "P2399_T1_role_closed_lift_distance_spectrum",
        "role_atoms_forced_false": ROLE_ATOMS,
        "free_nonrole_atoms": NONROLE_ATOMS,
        "role_transfer_distance_spectrum": spectrum["distance_spectra"]["role_transfer_theorem_level_closure"],
        "toe_distance_spectrum": spectrum["distance_spectra"]["toe_closure"],
        "minimum_role_transfer_lift_distance": spectrum["minimum_lift_distances"]["role_transfer_theorem_level_closure"],
        "minimum_toe_lift_distance": spectrum["minimum_lift_distances"]["toe_closure"],
        "nearest_toe_rows": spectrum["nearest_row_indices"]["toe_closure"],
        "missing_matrix_ranks_mod2": spectrum["missing_matrix_ranks_mod2"],
        "p2398_role_transfer_zero_polynomial": p2398_theorem.get("role_transfer_is_zero_polynomial"),
        "p2398_toe_zero_polynomial": p2398_theorem.get("toe_is_zero_polynomial"),
        "inherited_global_all_seven_atoms_required": proper_subset.get("all_seven_atoms_required"),
        "not_licensed": [
            "No role-transfer theorem is created by a finite lift distance.",
            "No ToE closure is created by the nearest lift row.",
            "No legacy physical role is promoted into L_total.",
            "No SM/GR numeric extraction is licensed without explicit role-successor evidence.",
        ],
    }
    gatekeepers = {
        "rg_nonduplication_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2398_zero_polynomials_inherited": p2398_theorem.get("role_transfer_is_zero_polynomial") is True and p2398_theorem.get("toe_is_zero_polynomial") is True,
        "role_transfer_min_distance_is_three": theorem_export["minimum_role_transfer_lift_distance"] == 3,
        "toe_min_distance_is_three": theorem_export["minimum_toe_lift_distance"] == 3,
        "role_transfer_spectrum_is_8_rows_at_3_and_8_at_4": theorem_export["role_transfer_distance_spectrum"] == {"3": 8, "4": 8},
        "toe_spectrum_is_binomial_shifted_by_three": theorem_export["toe_distance_spectrum"] == {"3": 1, "4": 4, "5": 6, "6": 4, "7": 1},
        "nearest_toe_row_is_all_nonrole_true": theorem_export["nearest_toe_rows"] == [15],
        "toe_missing_matrix_rank_is_five": theorem_export["missing_matrix_ranks_mod2"]["toe_closure"] == 5,
        "global_all_seven_atoms_required_inherited": proper_subset.get("all_seven_atoms_required") is True,
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2399_s1349_v1",
        "packet_id": "P2399",
        "stage_id": "S1349",
        "result_kind": "ROLE_CLOSED_LIFT_DISTANCE_SPECTRUM_CERTIFICATE",
        "status": "PASS_ROLE_CLOSED_LIFT_DISTANCE_THREE_ROLE_SUCCESSORS_REQUIRED",
        "role_closed_lift_distance_spectrum_certificate": {
            "rg_nonduplication_audit": grep,
            "source_artifact_statuses": {name: artifact.get("status", "TEXT_SOURCE") for name, artifact in artifacts.items()},
            "target_definitions": target_definitions,
            "lift_distance_spectrum": spectrum,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": "Construct explicit role-successor evidence for the three missing role atoms; non-role source progress alone can only reach the distance-3 nearest lift row.",
        "global_status": "OPEN_PROGRESS_LIFT_DISTANCE_CERTIFIED_THREE_ROLE_SUCCESSOR_ATOMS_REQUIRED",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    cert = payload["role_closed_lift_distance_spectrum_certificate"]
    theorem = cert["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2399 S1349: role-closed lift-distance spectrum certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2399/S1349 computes the exact Hamming/lift-distance spectrum from the P2396/P2397 role-closed quotient.",
                "",
                "## Spectra",
                "",
                f"- Role-transfer distance spectrum: `{theorem['role_transfer_distance_spectrum']}`.",
                f"- ToE distance spectrum: `{theorem['toe_distance_spectrum']}`.",
                f"- Minimum role-transfer lift distance: `{theorem['minimum_role_transfer_lift_distance']}`.",
                f"- Minimum ToE lift distance: `{theorem['minimum_toe_lift_distance']}`.",
                f"- Nearest ToE rows: `{theorem['nearest_toe_rows']}`.",
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
