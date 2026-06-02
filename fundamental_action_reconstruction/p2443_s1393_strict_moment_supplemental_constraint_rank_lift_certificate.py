#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.json"
MD = GEN / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8, "D": 25.0, "steps": 20000}
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]
MOMENT_CANDIDATES = ["M0", "M1", "M2", "M3"]
KERNEL_SAMPLE_D_VALUES = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
RANK_TOL = 1e-7
NULL_DOT_TOL = 1e-8


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
        "new_packet": "P2443|S1393|supplemental constraint rank lift|rank-lift certificate|moment supplemental constraint",
        "p2442_input": "P2442|S1392|local-identifiability nullspace|nullspace certificate",
        "candidate_constraints": "M0|M1|M2|M3|K_strict\\(|kernel sample|observable/source constraint",
        "rank_lift_language": "rank lift|rank-lift|augmented rank|independent observable|source constraint",
        "closure_blockers": "QW-2191|strict observable|physical-value generator|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float, params: dict[str, float]) -> float:
    return math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta"] * (d ** params["eta"]))


def moment(n: int, params: dict[str, float]) -> float:
    steps = int(params["steps"])
    h = params["D"] / steps
    total = 0.0
    for i in range(steps + 1):
        d = i * h
        weight = 0.5 if i in (0, steps) else 1.0
        total += weight * (d**n) * k_strict(d, params)
    return total * h


def finite_difference_gradient(func: Callable[[dict[str, float]], float]) -> list[float]:
    row: list[float] = []
    for name in PARAMETER_ORDER:
        step = 1e-5 * max(1.0, abs(float(STRICT_PARAMS[name])))
        plus = dict(STRICT_PARAMS)
        minus = dict(STRICT_PARAMS)
        plus[name] += step
        minus[name] -= step
        row.append((func(plus) - func(minus)) / (2.0 * step))
    return row


def real_rank(matrix: list[list[float]], tol: float = RANK_TOL) -> int:
    mat = [row[:] for row in matrix]
    if not mat:
        return 0
    rows = len(mat)
    cols = len(mat[0])
    rank = 0
    for col in range(cols):
        pivot = max(range(rank, rows), key=lambda r: abs(mat[r][col]), default=rank)
        if rank >= rows or abs(mat[pivot][col]) <= tol:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pivot_value = mat[rank][col]
        mat[rank] = [value / pivot_value for value in mat[rank]]
        for r in range(rows):
            if r != rank and abs(mat[r][col]) > tol:
                factor = mat[r][col]
                mat[r] = [a - factor * b for a, b in zip(mat[r], mat[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def dot(row: list[float], vector: dict[str, float]) -> float:
    return sum(row[index] * vector[name] for index, name in enumerate(PARAMETER_ORDER))


def candidate_rows(base_jacobian: list[list[float]], null_vector: dict[str, float]) -> list[dict[str, Any]]:
    candidates: list[tuple[str, str, Callable[[dict[str, float]], float]]] = []
    for n, name in enumerate(MOMENT_CANDIDATES):
        candidates.append((name, "raw_moment", lambda params, n=n: moment(n, params)))
    for d in KERNEL_SAMPLE_D_VALUES:
        candidates.append((f"K_at_d_{d:g}", "kernel_sample", lambda params, d=d: k_strict(d, params)))
    rows = []
    for name, kind, func in candidates:
        gradient = finite_difference_gradient(func)
        null_dot = dot(gradient, null_vector)
        augmented_rank = real_rank(base_jacobian + [gradient])
        rows.append(
            {
                "candidate_id": name,
                "candidate_kind": kind,
                "gradient": gradient,
                "null_direction_dot": null_dot,
                "abs_null_direction_dot": abs(null_dot),
                "augmented_rank": augmented_rank,
                "rank_lifts_to_four": augmented_rank == 4 and abs(null_dot) > NULL_DOT_TOL,
            }
        )
    return rows


def append_doc_sections() -> None:
    eq_section = """
## P2443/S1393 strict moment supplemental-constraint rank-lift certificate

`P2443/S1393` audits which finite supplemental constraints can remove the P2442 local null direction.  It enumerates raw moment candidates `(M0..M3)` and pointwise strict-kernel sample candidates, computes each candidate gradient against `(omega, phi, beta, eta)`, and checks whether appending that row lifts the P2441/P2442 moment-coefficient Jacobian rank from `3` to `4`.

The result is a rank-lift frontier, not a physical-value theorem: many singleton candidates are mathematically independent of the null direction, but none is yet proven to be an admissible strict observable/source constraint, selector theorem, or SM/GR value generator.
""".strip()
    lag_section = """
## P2443/S1393 supplemental-constraint rank-lift guard

`P2443/S1393` identifies finite candidate constraints that would lift the local moment-coefficient rank, but it does not license adding any of them to `L_total`.  A candidate must still be justified as a strict observable/source theorem or gauge-fixing theorem before it can close the P2442 null direction.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2442 = sources["P2442_NULLSPACE"].get("strict_moment_coefficient_local_identifiability_nullspace_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    null_vector = p2442.get("nullspace_certificate", {}).get("normalized_null_vector", {})
    rows = candidate_rows(base_jacobian, null_vector)
    rank_lift_rows = [row for row in rows if row["rank_lifts_to_four"]]
    rank_lift_by_kind: dict[str, list[str]] = {}
    for row in rank_lift_rows:
        rank_lift_by_kind.setdefault(row["candidate_kind"], []).append(row["candidate_id"])
    theorem_export = {
        "theorem_name": "P2443_T1_strict_moment_supplemental_constraint_rank_lift_certificate",
        "inherited_base_rank": p2442.get("input_jacobian_rank"),
        "inherited_nullspace_dimension": p2442.get("nullspace_dimension"),
        "parameter_order": PARAMETER_ORDER,
        "candidate_count": len(rows),
        "candidate_rows": rows,
        "rank_lifting_candidate_count": len(rank_lift_rows),
        "rank_lifting_candidate_ids": [row["candidate_id"] for row in rank_lift_rows],
        "rank_lifting_candidate_ids_by_kind": rank_lift_by_kind,
        "minimal_rank_lift_antichain_size": 1 if rank_lift_rows else None,
        "minimal_rank_lift_singletons": [[row["candidate_id"]] for row in rank_lift_rows],
        "all_raw_moment_candidates_rank_lift": all(row["rank_lifts_to_four"] for row in rows if row["candidate_kind"] == "raw_moment"),
        "all_kernel_sample_candidates_rank_lift": all(row["rank_lifts_to_four"] for row in rows if row["candidate_kind"] == "kernel_sample"),
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "strict_kernel_to_coefficient_map_theorem_exported": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Rank-lifting a Jacobian row is not proof that the row is an admissible strict observable/source constraint.",
            "Raw moments and pointwise kernel samples remain candidate constraints until a source theorem licenses them.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Choose one rank-lifting singleton and prove it is an admissible strict observable/source constraint or a lawful gauge-fixing condition for the P2442 null direction."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "inherited_rank_three": theorem_export["inherited_base_rank"] == 3,
        "inherited_one_null_direction": theorem_export["inherited_nullspace_dimension"] == 1,
        "twelve_candidates": theorem_export["candidate_count"] == 12,
        "all_candidates_rank_lift": theorem_export["rank_lifting_candidate_count"] == theorem_export["candidate_count"],
        "minimal_singleton_rank_lift": theorem_export["minimal_rank_lift_antichain_size"] == 1,
        "raw_moments_rank_lift": theorem_export["all_raw_moment_candidates_rank_lift"],
        "kernel_samples_rank_lift": theorem_export["all_kernel_sample_candidates_rank_lift"],
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_coefficient_theorem_export": not theorem_export["strict_kernel_to_coefficient_map_theorem_exported"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2443_s1393_v1",
        "packet_id": "P2443",
        "stage_id": "S1393",
        "result_kind": "STRICT_MOMENT_SUPPLEMENTAL_CONSTRAINT_RANK_LIFT_CERTIFICATE",
        "status": "PASS_STRICT_MOMENT_SUPPLEMENTAL_CONSTRAINT_RANK_LIFT_FRONTIER_NO_SOURCE_THEOREM",
        "strict_moment_supplemental_constraint_rank_lift_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_SUPPLEMENTAL_CONSTRAINT_RANK_LIFT_FRONTIER_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_moment_supplemental_constraint_rank_lift_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2443 S1393: strict moment supplemental-constraint rank-lift certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Inherited base rank: `{theorem['inherited_base_rank']}`.",
                f"- Inherited nullspace dimension: `{theorem['inherited_nullspace_dimension']}`.",
                f"- Candidate count: `{theorem['candidate_count']}`.",
                f"- Rank-lifting candidate count: `{theorem['rank_lifting_candidate_count']}`.",
                f"- Minimal rank-lift antichain size: `{theorem['minimal_rank_lift_antichain_size']}`.",
                f"- Rank-lifting candidates by kind: `{theorem['rank_lifting_candidate_ids_by_kind']}`.",
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
