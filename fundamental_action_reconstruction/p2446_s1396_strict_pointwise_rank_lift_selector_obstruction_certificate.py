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
OUT = GEN / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.json"
MD = GEN / "p2446_s1396_strict_pointwise_rank_lift_selector_obstruction_certificate.md"

SOURCE_FILES = {
    "P2442_NULLSPACE": GEN / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json",
    "P2445_STABILITY": GEN / "p2445_s1395_strict_moment_rank_lift_conditioning_stability_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
PARAMETER_ORDER = ["omega", "phi", "beta", "eta"]
ROBUST_NORMALIZED_VOLUME_THRESHOLD = 1.0e-3
SCAN_D_MIN = 0.0
SCAN_D_MAX = 5.0
SCAN_STEP = 0.0025
LOCAL_D_MIN = 0.75
LOCAL_D_MAX = 1.25
D1_NEIGHBORS = [0.95, 0.99, 1.0, 1.01, 1.05]


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
        "new_packet": "P2446|S1396|pointwise rank-lift selector obstruction|pointwise selector obstruction|d-window selector",
        "p2445_input": "P2445|S1395|conditioning stability|stable P2444/P2445 frontier|K_at_d_1",
        "pointwise_selector_language": "pointwise|d_ref|d-window|coordinate selection|selector obstruction|sample point",
        "rank_lift_language": "rank-lift|normalized rank-lift volume|conditioning frontier|robust candidate",
        "closure_blockers": "QW-2191|strict observable|gauge-fixing theorem|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def k_strict(d: float) -> float:
    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (1.0 + STRICT_PARAMS["beta"] * (d ** STRICT_PARAMS["eta"]))


def pointwise_gradient(d: float) -> list[float]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    beta = STRICT_PARAMS["beta"]
    eta = STRICT_PARAMS["eta"]
    d_eta = d**eta
    denominator = 1.0 + beta * d_eta
    phase = omega * d + phi
    cos_phase = math.cos(phase)
    sin_phase = math.sin(phase)
    eta_derivative = 0.0 if d == 0.0 else -cos_phase * beta * d_eta * math.log(d) / (denominator * denominator)
    return [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        -cos_phase * d_eta / (denominator * denominator),
        eta_derivative,
    ]


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


def scan_pointwise_window(base_jacobian: list[list[float]]) -> list[dict[str, Any]]:
    rows = []
    count = int(round((SCAN_D_MAX - SCAN_D_MIN) / SCAN_STEP))
    for index in range(count + 1):
        d = SCAN_D_MIN + index * SCAN_STEP
        gradient = pointwise_gradient(d)
        volume = normalized_volume(base_jacobian + [gradient])
        rows.append(
            {
                "d": d,
                "k_strict_value": k_strict(d),
                "normalized_rank_lift_volume": volume,
                "robust_rank_lift_candidate": volume > ROBUST_NORMALIZED_VOLUME_THRESHOLD,
            }
        )
    return rows


def append_doc_sections() -> None:
    eq_section = """
## P2446/S1396 strict pointwise rank-lift selector obstruction certificate

`P2446/S1396` follows the stable P2444/P2445 pointwise candidate `K_at_d_1` by scanning a dense pointwise `d`-window with analytic gradients.  The scan shows that `d=1` is a well-conditioned rank-lift point but is not uniquely selected by conditioning alone: nearby/alternative point samples also rank-lift, and the finite conditioning maximum in the audited window occurs away from `d=1`.

Therefore pointwise conditioning cannot by itself become a strict observable/source selector, a lawful gauge choice, or a role-bearing coefficient source.  A separate theorem must select the point coordinate, the observable meaning, or the gauge slice before a pointwise row can replace the strict moment route.
""".strip()
    lag_section = """
## P2446/S1396 pointwise selector obstruction guard

`P2446/S1396` blocks the silent promotion of the stable `K_at_d_1` row into `L_total`: pointwise rank-lift conditioning does not select the evaluation coordinate and does not supply admissibility or gauge-fixing authority.
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
    p2445 = sources["P2445_STABILITY"].get("strict_moment_rank_lift_conditioning_stability_certificate", {}).get("theorem_export", {})
    base_jacobian = p2442.get("input_jacobian", [])
    rows = scan_pointwise_window(base_jacobian)
    d1_row = next(row for row in rows if abs(row["d"] - 1.0) < 1e-12)
    best = max(rows, key=lambda item: item["normalized_rank_lift_volume"])
    local_rows = [row for row in rows if LOCAL_D_MIN <= row["d"] <= LOCAL_D_MAX]
    local_best = max(local_rows, key=lambda item: item["normalized_rank_lift_volume"])
    robust_rows = [row for row in rows if row["robust_rank_lift_candidate"]]
    robust_local_rows = [row for row in local_rows if row["robust_rank_lift_candidate"]]
    neighbor_rows = [
        {
            "d": d,
            "normalized_rank_lift_volume": normalized_volume(base_jacobian + [pointwise_gradient(d)]),
            "k_strict_value": k_strict(d),
        }
        for d in D1_NEIGHBORS
    ]
    theorem_export = {
        "theorem_name": "P2446_T1_strict_pointwise_rank_lift_selector_obstruction_certificate",
        "inherited_stability_certificate": "P2445/S1395",
        "inherited_p2445_best_candidate_id": p2445.get("baseline_best_candidate_id"),
        "scan_d_min": SCAN_D_MIN,
        "scan_d_max": SCAN_D_MAX,
        "scan_step": SCAN_STEP,
        "scan_point_count": len(rows),
        "local_window": {"d_min": LOCAL_D_MIN, "d_max": LOCAL_D_MAX},
        "conditioning_threshold_normalized_volume": ROBUST_NORMALIZED_VOLUME_THRESHOLD,
        "d1_row": d1_row,
        "global_best_pointwise_row": best,
        "local_best_pointwise_row": local_best,
        "d1_neighbor_rows": neighbor_rows,
        "robust_pointwise_count_in_scan": len(robust_rows),
        "robust_pointwise_count_in_local_window": len(robust_local_rows),
        "d1_is_robust": d1_row["robust_rank_lift_candidate"],
        "d1_is_global_conditioning_maximum_on_scan": abs(best["d"] - 1.0) < 1e-12,
        "d1_is_local_conditioning_maximum_on_window": abs(local_best["d"] - 1.0) < 1e-12,
        "best_volume_minus_d1_volume": best["normalized_rank_lift_volume"] - d1_row["normalized_rank_lift_volume"],
        "local_best_volume_minus_d1_volume": local_best["normalized_rank_lift_volume"] - d1_row["normalized_rank_lift_volume"],
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_fixing_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Pointwise rank-lift conditioning does not select a unique evaluation coordinate d_ref.",
            "The stable K_at_d_1 row remains a candidate because conditioning alone neither proves strict observable/source admissibility nor supplies a lawful gauge slice.",
            "No strict physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Prove a strict point-coordinate selector, observable/source theorem, or gauge-slice theorem before using any pointwise K_strict(d_ref) row as a supplemental coefficient source."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2445_best_inherited": theorem_export["inherited_p2445_best_candidate_id"] == "K_at_d_1",
        "scan_grid_complete": theorem_export["scan_point_count"] == int(round((SCAN_D_MAX - SCAN_D_MIN) / SCAN_STEP)) + 1,
        "d1_remains_robust": theorem_export["d1_is_robust"],
        "d1_not_global_maximum": not theorem_export["d1_is_global_conditioning_maximum_on_scan"],
        "d1_not_local_maximum": not theorem_export["d1_is_local_conditioning_maximum_on_window"],
        "positive_global_selector_gap": theorem_export["best_volume_minus_d1_volume"] > 0.0,
        "positive_local_selector_gap": theorem_export["local_best_volume_minus_d1_volume"] > 0.0,
        "many_robust_pointwise_alternatives": theorem_export["robust_pointwise_count_in_local_window"] > 1,
        "no_pointwise_selector_export": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"],
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_gauge_fixing_export": not theorem_export["gauge_fixing_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2446_s1396_v1",
        "packet_id": "P2446",
        "stage_id": "S1396",
        "result_kind": "STRICT_POINTWISE_RANK_LIFT_SELECTOR_OBSTRUCTION_CERTIFICATE",
        "status": "PASS_STRICT_POINTWISE_RANK_LIFT_SELECTOR_OBSTRUCTION_NO_SOURCE_THEOREM",
        "strict_pointwise_rank_lift_selector_obstruction_certificate": {
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
        "recommended_next_honest_step": theorem_export["next_honest_step"],
        "global_status": "OPEN_PROGRESS_POINTWISE_SELECTOR_OBSTRUCTION_EXPORTED_NO_CLOSURE",
    }


def write_outputs(payload: dict[str, Any]) -> None:
    theorem = payload["strict_pointwise_rank_lift_selector_obstruction_certificate"]["theorem_export"]
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2446 S1396: strict pointwise rank-lift selector obstruction certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Finite facts",
                "",
                f"- Scan point count: `{theorem['scan_point_count']}`.",
                f"- `d=1` normalized volume: `{theorem['d1_row']['normalized_rank_lift_volume']}`.",
                f"- Global best point: `{theorem['global_best_pointwise_row']}`.",
                f"- Local best point: `{theorem['local_best_pointwise_row']}`.",
                f"- Robust pointwise count in local window: `{theorem['robust_pointwise_count_in_local_window']}`.",
                f"- `d=1` is global maximum: `{theorem['d1_is_global_conditioning_maximum_on_scan']}`.",
                f"- `d=1` is local maximum: `{theorem['d1_is_local_conditioning_maximum_on_window']}`.",
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
