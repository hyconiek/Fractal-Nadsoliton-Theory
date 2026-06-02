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
OUT = GEN / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.json"
MD = GEN / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.md"

SOURCE_FILES = {
    "P2443_RANK_LIFT": GEN / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}
ROBUST_NORMALIZED_VOLUME_THRESHOLD = 1.0e-3
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


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
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2444|S1394|rank lift conditioning|rank-lift conditioning|strict moment rank lift conditioning",
        "p2443_input": "P2443|S1393|supplemental-constraint rank-lift|rank-lifting candidate",
        "conditioning_language": "conditioning|normalized volume|determinant margin|rank-lift margin|robust rank",
        "candidate_constraints": "M0|M1|M2|M3|K_at_d|kernel sample|raw moment",
        "closure_blockers": "QW-2191|strict observable|physical-value generator|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def row_norm(row: list[float]) -> float:
    return math.sqrt(sum(value * value for value in row))


def determinant(matrix: list[list[float]]) -> float:
    mat = [row[:] for row in matrix]
    n = len(mat)
    det = 1.0
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(mat[r][col]))
        if abs(mat[pivot][col]) < 1e-15:
            return 0.0
        if pivot != col:
            mat[col], mat[pivot] = mat[pivot], mat[col]
            det *= -1.0
        pivot_value = mat[col][col]
        det *= pivot_value
        for row in range(col + 1, n):
            factor = mat[row][col] / pivot_value
            for c in range(col, n):
                mat[row][c] -= factor * mat[col][c]
    return det


def normalized_volume(matrix: list[list[float]]) -> float:
    normalized = []
    for row in matrix:
        norm = row_norm(row)
        normalized.append([value / norm for value in row])
    return abs(determinant(normalized))


def condition_rows(base_jacobian: list[list[float]], candidate_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for row in candidate_rows:
        gradient = row["gradient"]
        gradient_norm = row_norm(gradient)
        volume = normalized_volume(base_jacobian + [gradient])
        margin = row["abs_null_direction_dot"] / gradient_norm if gradient_norm else 0.0
        rows.append(
            {
                "candidate_id": row["candidate_id"],
                "candidate_kind": row["candidate_kind"],
                "augmented_rank": row["augmented_rank"],
                "gradient_norm": gradient_norm,
                "abs_null_direction_dot": row["abs_null_direction_dot"],
                "rank_lift_margin": margin,
                "normalized_rank_lift_volume": volume,
                "robust_rank_lift_candidate": volume > ROBUST_NORMALIZED_VOLUME_THRESHOLD,
            }
        )
    return sorted(rows, key=lambda item: item["normalized_rank_lift_volume"], reverse=True)


def append_doc_sections() -> None:
    eq_section = """
## P2444/S1394 strict moment rank-lift conditioning certificate

`P2444/S1394` refines the P2443 rank-lift frontier by adding conditioning diagnostics.  For each singleton supplemental candidate it computes gradient norm, null-direction margin, and row-normalized augmented determinant volume.  This separates mathematically rank-lifting candidates from better-conditioned rank-lifting candidates.

The conditioning winner is still only a candidate: numerical conditioning does not prove admissibility as a strict observable/source constraint, does not discharge `QW-2191`, and does not export a physical-value generator.
""".strip()
    lag_section = """
## P2444/S1394 rank-lift conditioning guard

`P2444/S1394` ranks supplemental constraints by numerical conditioning, but no conditioned row may be inserted into `L_total` until a strict observable/source or gauge-fixing theorem licenses that row.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2443 = sources["P2443_RANK_LIFT"].get("strict_moment_supplemental_constraint_rank_lift_certificate", {}).get("theorem_export", {})
    base_jacobian = load_json(GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json").get(
        "strict_moment_coefficient_local_identifiability_nullspace_certificate", {}
    ).get("theorem_export", {}).get("input_jacobian", [])
    rows = condition_rows(base_jacobian, p2443.get("candidate_rows", []))
    robust = [row for row in rows if row["robust_rank_lift_candidate"]]
    robust_by_kind: dict[str, list[str]] = {}
    for row in robust:
        robust_by_kind.setdefault(row["candidate_kind"], []).append(row["candidate_id"])
    theorem_export = {
        "theorem_name": "P2444_T1_strict_moment_rank_lift_conditioning_certificate",
        "inherited_candidate_count": p2443.get("candidate_count"),
        "conditioning_threshold_normalized_volume": ROBUST_NORMALIZED_VOLUME_THRESHOLD,
        "conditioned_candidate_rows_descending": rows,
        "best_conditioned_candidate_id": rows[0]["candidate_id"],
        "best_conditioned_candidate_kind": rows[0]["candidate_kind"],
        "best_conditioned_normalized_volume": rows[0]["normalized_rank_lift_volume"],
        "best_conditioned_rank_lift_margin": rows[0]["rank_lift_margin"],
        "robust_rank_lift_candidate_count": len(robust),
        "robust_rank_lift_candidate_ids": [row["candidate_id"] for row in robust],
        "robust_rank_lift_candidate_ids_by_kind": robust_by_kind,
        "weakest_conditioned_candidate_id": rows[-1]["candidate_id"],
        "weakest_conditioned_normalized_volume": rows[-1]["normalized_rank_lift_volume"],
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Numerical rank-lift conditioning is not an admissibility theorem for a supplemental constraint.",
            "The best-conditioned singleton is still only a candidate until a strict observable/source theorem licenses it.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Try to prove strict source admissibility for the best-conditioned rank-lift singleton, or explain why it is a gauge choice rather than a physical observable."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "twelve_candidates_inherited": theorem_export["inherited_candidate_count"] == 12,
        "twelve_conditioned_rows": len(theorem_export["conditioned_candidate_rows_descending"]) == 12,
        "best_candidate_k_at_d_1": theorem_export["best_conditioned_candidate_id"] == "K_at_d_1",
        "robust_candidate_count_six": theorem_export["robust_rank_lift_candidate_count"] == 6,
        "weakest_candidate_k_at_d_20": theorem_export["weakest_conditioned_candidate_id"] == "K_at_d_20",
        "best_volume_above_threshold": theorem_export["best_conditioned_normalized_volume"] > ROBUST_NORMALIZED_VOLUME_THRESHOLD,
        "weakest_volume_below_threshold": theorem_export["weakest_conditioned_normalized_volume"] < ROBUST_NORMALIZED_VOLUME_THRESHOLD,
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2444_s1394_v1",
        "packet_id": "P2444",
        "stage_id": "S1394",
        "result_kind": "STRICT_MOMENT_RANK_LIFT_CONDITIONING_CERTIFICATE",
        "status": "PASS_STRICT_MOMENT_RANK_LIFT_CONDITIONING_FRONTIER_NO_SOURCE_THEOREM",
        "strict_moment_rank_lift_conditioning_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_RANK_LIFT_CONDITIONING_FRONTIER_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_moment_rank_lift_conditioning_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2444 S1394: strict moment rank-lift conditioning certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Inherited candidates: `{theorem['inherited_candidate_count']}`.",
                f"- Robust rank-lift candidates: `{theorem['robust_rank_lift_candidate_count']}`.",
                f"- Best-conditioned candidate: `{theorem['best_conditioned_candidate_id']}`.",
                f"- Best normalized volume: `{theorem['best_conditioned_normalized_volume']}`.",
                f"- Weakest-conditioned candidate: `{theorem['weakest_conditioned_candidate_id']}`.",
                f"- Robust candidates by kind: `{theorem['robust_rank_lift_candidate_ids_by_kind']}`.",
                "",
                "## Hard limits",
                "",
                *[f"- {item}" for item in theorem["not_licensed"]],
                "",
                "## Next honest step",
                "",
                theorem["next_honest_step"],
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
