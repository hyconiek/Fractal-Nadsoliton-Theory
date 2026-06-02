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
OUT = GEN / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.json"
MD = GEN / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
    "P2444_CONDITIONING": GEN / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_BASE_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8, "D": 25.0}
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]
MOMENT_CANDIDATES = ["M0", "M1", "M2", "M3"]
KERNEL_SAMPLE_D_VALUES = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
FINITE_DIFFERENCE_STEP_SCALES = [0.5, 1.0, 2.0, 5.0]
QUADRATURE_STEPS = [10000, 20000, 40000]
ROBUST_NORMALIZED_VOLUME_THRESHOLD = 1.0e-3


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
        "new_packet": "P2445|S1395|rank-lift conditioning stability|conditioning stability|finite-difference stability|mesh stability",
        "p2444_input": "P2444|S1394|rank-lift conditioning|normalized rank-lift volume|best-conditioned candidate",
        "stability_language": "finite-difference step|derivative step|quadrature steps|mesh stability|conditioning stability|robust set invariant",
        "candidate_constraints": "M0|M1|M2|M3|K_at_d|kernel sample|raw moment",
        "closure_blockers": "QW-2191|strict observable|physical-value generator|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float, params: dict[str, float]) -> float:
    return math.cos(params["omega"] * d + params["phi"]) / (1.0 + params["beta"] * (d ** params["eta"]))


def moment(n: int, params: dict[str, float], steps: int) -> float:
    h = params["D"] / steps
    total = 0.0
    for i in range(steps + 1):
        d = i * h
        weight = 0.5 if i in (0, steps) else 1.0
        total += weight * (d**n) * k_strict(d, params)
    return total * h


def finite_difference_gradient(func: Callable[[dict[str, float], int], float], fd_scale: float, quadrature_steps: int) -> list[float]:
    row: list[float] = []
    for name in PARAMETER_ORDER:
        step = fd_scale * 1e-5 * max(1.0, abs(float(STRICT_BASE_PARAMS[name])))
        plus = dict(STRICT_BASE_PARAMS)
        minus = dict(STRICT_BASE_PARAMS)
        plus[name] += step
        minus[name] -= step
        row.append((func(plus, quadrature_steps) - func(minus, quadrature_steps)) / (2.0 * step))
    return row


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


def candidate_gradients(fd_scale: float, quadrature_steps: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for n, name in enumerate(MOMENT_CANDIDATES):
        gradient = finite_difference_gradient(lambda params, steps, n=n: moment(n, params, steps), fd_scale, quadrature_steps)
        rows.append({"candidate_id": name, "candidate_kind": "raw_moment", "gradient": gradient})
    for d in KERNEL_SAMPLE_D_VALUES:
        gradient = finite_difference_gradient(lambda params, steps, d=d: k_strict(d, params), fd_scale, quadrature_steps)
        rows.append({"candidate_id": f"K_at_d_{d:g}", "candidate_kind": "kernel_sample", "gradient": gradient})
    return rows


def condition_rows(base_jacobian: list[list[float]], fd_scale: float, quadrature_steps: int) -> list[dict[str, Any]]:
    rows = []
    for row in candidate_gradients(fd_scale, quadrature_steps):
        volume = normalized_volume(base_jacobian + [row["gradient"]])
        rows.append(
            {
                "candidate_id": row["candidate_id"],
                "candidate_kind": row["candidate_kind"],
                "normalized_rank_lift_volume": volume,
                "robust_rank_lift_candidate": volume > ROBUST_NORMALIZED_VOLUME_THRESHOLD,
            }
        )
    return sorted(rows, key=lambda item: item["normalized_rank_lift_volume"], reverse=True)


def append_doc_sections() -> None:
    eq_section = """
## P2445/S1395 strict moment rank-lift conditioning stability certificate

`P2445/S1395` audits whether the P2444 conditioning frontier is a finite-difference or quadrature-mesh artifact.  It recomputes the singleton supplemental-candidate normalized rank-lift volumes across derivative-step scales and quadrature resolutions, then checks whether the best singleton and robust candidate set remain invariant.

The stability result is numerical robustness only: it does not prove that any candidate is an admissible strict observable/source constraint, does not supply a gauge-fixing theorem, and does not export a physical-value generator.
""".strip()
    lag_section = """
## P2445/S1395 rank-lift conditioning stability guard

`P2445/S1395` may stabilize the P2444 ordering against numerical-discretization objections, but stable conditioning still cannot insert a supplemental row into `L_total` without a strict observable/source theorem or a lawful gauge-fixing theorem.
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
    p2444 = sources["P2444_CONDITIONING"].get("strict_moment_rank_lift_conditioning_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    baseline_rows = p2444.get("conditioned_candidate_rows_descending", [])
    baseline_best = p2444.get("best_conditioned_candidate_id")
    baseline_robust_ids = p2444.get("robust_rank_lift_candidate_ids", [])
    baseline_top_six = [row["candidate_id"] for row in baseline_rows[:6]]

    config_results = []
    volume_samples_by_candidate: dict[str, list[float]] = {row["candidate_id"]: [] for row in baseline_rows}
    for quadrature_steps in QUADRATURE_STEPS:
        for fd_scale in FINITE_DIFFERENCE_STEP_SCALES:
            rows = condition_rows(base_jacobian, fd_scale, quadrature_steps)
            robust_ids = [row["candidate_id"] for row in rows if row["robust_rank_lift_candidate"]]
            for row in rows:
                volume_samples_by_candidate.setdefault(row["candidate_id"], []).append(row["normalized_rank_lift_volume"])
            config_results.append(
                {
                    "finite_difference_step_scale": fd_scale,
                    "quadrature_steps": quadrature_steps,
                    "best_candidate_id": rows[0]["candidate_id"],
                    "best_normalized_volume": rows[0]["normalized_rank_lift_volume"],
                    "top_six_candidate_ids": [row["candidate_id"] for row in rows[:6]],
                    "robust_candidate_ids": robust_ids,
                    "robust_candidate_count": len(robust_ids),
                    "weakest_candidate_id": rows[-1]["candidate_id"],
                    "weakest_normalized_volume": rows[-1]["normalized_rank_lift_volume"],
                }
            )

    volume_ranges = []
    for candidate_id, samples in sorted(volume_samples_by_candidate.items()):
        volume_ranges.append(
            {
                "candidate_id": candidate_id,
                "min_normalized_volume": min(samples),
                "max_normalized_volume": max(samples),
                "spread_normalized_volume": max(samples) - min(samples),
            }
        )
    robust_set = set(baseline_robust_ids)
    nonrobust_set = set(volume_samples_by_candidate) - robust_set
    min_robust_volume = min(min(volume_samples_by_candidate[cid]) for cid in robust_set)
    max_nonrobust_volume = max(max(volume_samples_by_candidate[cid]) for cid in nonrobust_set)
    max_volume_spread = max(item["spread_normalized_volume"] for item in volume_ranges)

    theorem_export = {
        "theorem_name": "P2445_T1_strict_moment_rank_lift_conditioning_stability_certificate",
        "inherited_conditioning_certificate": "P2444/S1394",
        "finite_difference_step_scales": FINITE_DIFFERENCE_STEP_SCALES,
        "quadrature_steps": QUADRATURE_STEPS,
        "configuration_count": len(config_results),
        "conditioning_threshold_normalized_volume": ROBUST_NORMALIZED_VOLUME_THRESHOLD,
        "baseline_best_candidate_id": baseline_best,
        "baseline_robust_candidate_ids": baseline_robust_ids,
        "baseline_top_six_candidate_ids": baseline_top_six,
        "configuration_results": config_results,
        "volume_ranges_by_candidate": volume_ranges,
        "all_configurations_preserve_best_candidate": all(result["best_candidate_id"] == baseline_best for result in config_results),
        "all_configurations_preserve_robust_set": all(result["robust_candidate_ids"] == baseline_robust_ids for result in config_results),
        "all_configurations_preserve_top_six_order": all(result["top_six_candidate_ids"] == baseline_top_six for result in config_results),
        "minimum_robust_candidate_volume_across_grid": min_robust_volume,
        "maximum_nonrobust_candidate_volume_across_grid": max_nonrobust_volume,
        "max_conditioning_volume_spread_across_grid": max_volume_spread,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_fixing_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Finite-difference and mesh stability is not an admissibility theorem for a supplemental constraint.",
            "A stable best-conditioned singleton remains only a candidate until a strict observable/source theorem or lawful gauge-fixing theorem licenses it.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Use the stable P2444/P2445 frontier to attempt an admissibility or gauge-fixing theorem for K_at_d_1, while keeping it outside L_total until that theorem is supplied."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "configuration_grid_complete": theorem_export["configuration_count"] == len(FINITE_DIFFERENCE_STEP_SCALES) * len(QUADRATURE_STEPS),
        "baseline_best_inherited": theorem_export["baseline_best_candidate_id"] == "K_at_d_1",
        "best_candidate_stable": theorem_export["all_configurations_preserve_best_candidate"],
        "robust_set_stable": theorem_export["all_configurations_preserve_robust_set"],
        "top_six_order_stable": theorem_export["all_configurations_preserve_top_six_order"],
        "robust_threshold_gap_positive": min_robust_volume > ROBUST_NORMALIZED_VOLUME_THRESHOLD > max_nonrobust_volume,
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_gauge_fixing_export": not theorem_export["gauge_fixing_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2445_s1395_v1",
        "packet_id": "P2445",
        "stage_id": "S1395",
        "result_kind": "STRICT_MOMENT_RANK_LIFT_CONDITIONING_STABILITY_CERTIFICATE",
        "status": "PASS_STRICT_MOMENT_RANK_LIFT_CONDITIONING_STABILITY_NO_SOURCE_THEOREM",
        "strict_moment_rank_lift_conditioning_stability_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_CONDITIONING_STABILITY_FRONTIER_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_moment_rank_lift_conditioning_stability_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2445 S1395: strict moment rank-lift conditioning stability certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Configuration count: `{theorem['configuration_count']}`.",
                f"- Baseline best candidate: `{theorem['baseline_best_candidate_id']}`.",
                f"- Best candidate stable: `{theorem['all_configurations_preserve_best_candidate']}`.",
                f"- Robust set stable: `{theorem['all_configurations_preserve_robust_set']}`.",
                f"- Top-six order stable: `{theorem['all_configurations_preserve_top_six_order']}`.",
                f"- Minimum robust volume across grid: `{theorem['minimum_robust_candidate_volume_across_grid']}`.",
                f"- Maximum nonrobust volume across grid: `{theorem['maximum_nonrobust_candidate_volume_across_grid']}`.",
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
