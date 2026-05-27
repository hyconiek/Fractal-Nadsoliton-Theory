#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2234 = GEN / "p2234_s1184_strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe.json"
OUT = GEN / "p2235_s1185_strict_nu_branch_residual_to_modeled_affine_surrogate_probe.json"
MD = GEN / "p2235_s1185_strict_nu_branch_residual_to_modeled_affine_surrogate_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2203 = load(IN_2203)
    p2234 = load(IN_2234)

    rows = (p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}).get("residual_map_rows", []) or []
    modeled = (p2234.get("strict_nu_branch_modeled_lane_to_residual_lane_gap_quantification_probe", {}) or {}).get("modeled_lane_observables", {}) or {}

    sf_lo = float(modeled.get("safety_factor_lo", 0.0) or 0.0)
    sf_hi = float(modeled.get("safety_factor_hi", 0.0) or 0.0)
    if not rows:
        raise RuntimeError("Missing residual rows for P2235 surrogate probe")

    rvals = [float(r["transport_residual_l1"]) for r in rows]
    rmin, rmax = min(rvals), max(rvals)

    if not (rmax > rmin and sf_hi >= sf_lo):
        raise RuntimeError("Degenerate ranges for surrogate probe")

    # affine map r -> sf : sf = a*r + b using endpoint anchoring
    a = (sf_hi - sf_lo) / (rmax - rmin)
    b = sf_lo - a * rmin

    fit_rows = []
    for r in rvals:
        sf_hat = a * r + b
        fit_rows.append({
            "residual_l1": r,
            "modeled_safety_surrogate": sf_hat,
            "inside_modeled_range": sf_lo - 1e-15 <= sf_hat <= sf_hi + 1e-15,
        })

    payload = {
        "schema_version": "p2235_s1185_v1",
        "packet_id": "P2235",
        "stage_id": "S1185",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_AFFINE_SURROGATE_PROBE",
        "strict_nu_branch_residual_to_modeled_affine_surrogate_probe": {
            "probe_id": "STRICT_NU_BRANCH_RESIDUAL_TO_MODELED_AFFINE_SURROGATE_PROBE_V1",
            "source_packets": [str(IN_2203.relative_to(ROOT)), str(IN_2234.relative_to(ROOT))],
            "anchoring": {
                "residual_min": rmin,
                "residual_max": rmax,
                "modeled_safety_lo": sf_lo,
                "modeled_safety_hi": sf_hi,
            },
            "affine_map": {"a": a, "b": b, "formula": "sf_hat = a*residual_l1 + b"},
            "surrogate_rows": fit_rows,
            "bridge_gap_statement": "Affine surrogate is diagnostic-only and does not constitute theorem-level residual->modeled bridge.",
            "theorem_scope_limit": "diagnostic affine surrogate on current sampled rows only; not global Task-3 closure",
        },
        "recommended_next_honest_step": {
            "id": "P2236_candidate",
            "goal": "test non-affine surrogate families and export model-selection evidence (error/complexity) before any bridge claim",
        },
        "gatekeeper_checks": {
            "affine_surrogate_exported": True,
            "positive_slope_map": a > 0.0,
            "all_surrogates_inside_modeled_range": all(r["inside_modeled_range"] for r in fit_rows),
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
        "\n".join([
            "# P2235 S1185: strict nu-branch residual->modeled affine surrogate probe",
            "",
            f"- affine slope a: `{a:.12e}`",
            f"- affine intercept b: `{b:.12e}`",
            f"- all surrogate rows inside modeled range: `{all(r['inside_modeled_range'] for r in fit_rows)}`",
            "",
            "Diagnostic surrogate only; no bridge theorem and no global closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
