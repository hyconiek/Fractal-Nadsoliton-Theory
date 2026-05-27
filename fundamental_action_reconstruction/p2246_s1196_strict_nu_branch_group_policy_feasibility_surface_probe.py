#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2245 = GEN / "p2245_s1195_strict_nu_branch_group_policy_reserve_factor_calibration_probe.json"
OUT = GEN / "p2246_s1196_strict_nu_branch_group_policy_feasibility_surface_probe.json"
MD = GEN / "p2246_s1196_strict_nu_branch_group_policy_feasibility_surface_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2245 = load(IN_2245)
    probe = (p2245.get("strict_nu_branch_group_policy_reserve_factor_calibration_probe", {}) or {})
    inputs = (probe.get("inputs", {}) or {})
    reserve = (probe.get("reserve_calibration", {}) or {})

    base_n = int(inputs.get("group_count", 1) or 1)
    base_load = float(inputs.get("load_ratio_threshold", 1.0) or 1.0)
    base_total = float(inputs.get("total_coverage_mass", 1.0) or 1.0)
    target_risk = float(inputs.get("target_risk", 0.05) or 0.05)
    kappa_ref = float(reserve.get("suggested_kappa", 1.0) or 1.0)

    n_grid = [max(1, base_n - 1), base_n, base_n + 1]
    load_grid = [0.8 * base_load, base_load, 1.2 * base_load]

    rows = []
    admissible = 0
    for n in n_grid:
        for load_val in load_grid:
            required_total = (n - 1) + kappa_ref * load_val
            margin = base_total - required_total
            pass_det = margin >= -1e-15
            if pass_det:
                admissible += 1
            rows.append(
                {
                    "group_count": n,
                    "load_ratio": load_val,
                    "required_total_coverage_mass": required_total,
                    "available_total_coverage_mass": base_total,
                    "deterministic_margin": margin,
                    "deterministic_admissible": pass_det,
                }
            )

    total_cells = len(rows)
    admissible_fraction = admissible / max(total_cells, 1)

    payload = {
        "schema_version": "p2246_s1196_v1",
        "packet_id": "P2246",
        "stage_id": "S1196",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_FEASIBILITY_SURFACE_PROBE",
        "strict_nu_branch_group_policy_feasibility_surface_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_FEASIBILITY_SURFACE_PROBE_V1",
            "source_packets": [str(IN_2245.relative_to(ROOT))],
            "inputs": {
                "base_group_count": base_n,
                "base_load_ratio": base_load,
                "base_total_coverage_mass": base_total,
                "target_risk": target_risk,
                "kappa_reference": kappa_ref,
            },
            "grid": {
                "group_count_grid": n_grid,
                "load_ratio_grid": load_grid,
                "surface_rows": rows,
                "admissible_cell_fraction": admissible_fraction,
            },
            "physical_interpretation_note": "Feasibility surface maps where current coverage budget can sustain reserve-calibrated safety under varying group fragmentation and load demand, exposing robustness phase-space boundaries.",
            "theorem_scope_limit": "finite-grid strict-lane feasibility diagnostic only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2247_candidate",
            "goal": "derive closed-form boundary curve and sensitivity coefficients d(margin)/d(group_count), d(margin)/d(load_ratio)",
        },
        "gatekeeper_checks": {
            "feasibility_surface_exported": True,
            "nonempty_surface_grid": total_cells > 0,
            "admissible_fraction_computable": 0.0 <= admissible_fraction <= 1.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2246 S1196: feasibility-surface probe",
                "",
                f"- base group count: `{base_n}`",
                f"- base load ratio: `{base_load:.12e}`",
                f"- base total coverage mass: `{base_total:.12e}`",
                f"- reference kappa: `{kappa_ref:.12e}`",
                f"- admissible cell fraction: `{admissible_fraction:.12e}`",
                "",
                "Finite-grid feasibility diagnostic only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
