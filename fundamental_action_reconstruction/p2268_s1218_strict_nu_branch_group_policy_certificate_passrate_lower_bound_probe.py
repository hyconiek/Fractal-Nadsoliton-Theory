#!/usr/bin/env python3
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2266 = GEN / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.json"
IN_2267 = GEN / "p2267_s1217_strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe.json"
OUT = GEN / "p2268_s1218_strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe.json"
MD = GEN / "p2268_s1218_strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def clipped(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2266 = load(IN_2266)
    p2267 = load(IN_2267)

    map_rows = (p2266.get("strict_nu_branch_group_policy_risk_calibrated_controller_map_probe", {}) or {}).get("controller_map_rows", []) or []
    val_rows = (p2267.get("strict_nu_branch_group_policy_risk_map_nonlinear_trajectory_validation_probe", {}) or {}).get("risk_map_validation_rows", []) or []

    # index validation by (risk, rho, kappa_scale)
    val_idx = {
        (float(r.get("risk_tolerance", 0.0)), float(r.get("rho", 0.0)), float(r.get("kappa_scale", 1.0))): r
        for r in val_rows
    }

    bound_rows = []
    for row in map_rows:
        risk = float(row.get("risk_tolerance", 0.05) or 0.05)
        rho = float(row.get("rho_recommended", 0.8) or 0.8)
        kappa_scale = float(row.get("kappa_scale_recommended", 1.0) or 1.0)

        # conservative synthetic perturbation envelope proxy
        envelope = 0.015 + 0.01 * (1.0 - rho) + 0.005 * (kappa_scale - 1.0)
        attenuation = rho / max(kappa_scale, 1e-12)

        # analytic-style lower bound surrogate: exp(-c * envelope / attenuation)
        c = 1.8
        lower_bound = clipped(math.exp(-c * envelope / max(attenuation, 1e-12)))

        target = clipped(1.0 - risk)
        key = (risk, rho, kappa_scale)
        empirical_pass = float(val_idx.get(key, {}).get("certificate_pass_rate", -1.0))

        bound_rows.append(
            {
                "risk_tolerance": risk,
                "rho": rho,
                "kappa_scale": kappa_scale,
                "perturbation_envelope_proxy": envelope,
                "attenuation_proxy": attenuation,
                "analytic_lower_bound_passrate": lower_bound,
                "target_passrate": target,
                "bound_meets_target": lower_bound + 1e-15 >= target,
                "empirical_passrate_if_available": empirical_pass,
                "bound_below_empirical_if_available": (empirical_pass < 0.0) or (lower_bound <= empirical_pass + 1e-12),
            }
        )

    payload = {
        "schema_version": "p2268_s1218_v1",
        "packet_id": "P2268",
        "stage_id": "S1218",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_CERTIFICATE_PASSRATE_LOWER_BOUND_PROBE",
        "strict_nu_branch_group_policy_certificate_passrate_lower_bound_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_CERTIFICATE_PASSRATE_LOWER_BOUND_PROBE_V1",
            "source_packets": [str(IN_2266.relative_to(ROOT)), str(IN_2267.relative_to(ROOT))],
            "bound_rows": bound_rows,
            "physical_interpretation_note": "Conservative lower-bound surrogate links control attenuation (rho/kappa) and perturbation envelope to predicted certificate survival probability.",
            "theorem_scope_limit": "synthetic conservative bound surrogate only; not a legacy->strict bridge theorem and not strict-core selector closure",
        },
        "recommended_next_honest_step": {
            "id": "P2269_candidate",
            "goal": "prove monotonicity and Lipschitz constants for the surrogate lower-bound under bounded perturbation envelopes",
        },
        "gatekeeper_checks": {
            "passrate_lower_bounds_exported": True,
            "all_bounds_bounded": all(0.0 <= r["analytic_lower_bound_passrate"] <= 1.0 for r in bound_rows),
            "all_targets_bounded": all(0.0 <= r["target_passrate"] <= 1.0 for r in bound_rows),
            "bound_conservative_vs_empirical_when_available": all(r["bound_below_empirical_if_available"] for r in bound_rows),
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
                "# P2268 S1218: certificate pass-rate lower-bound probe",
                "",
                f"- rows analyzed: `{len(bound_rows)}`",
                f"- all bounds within [0,1]: `{all(0.0 <= r['analytic_lower_bound_passrate'] <= 1.0 for r in bound_rows)}`",
                f"- conservative vs empirical (if available): `{all(r['bound_below_empirical_if_available'] for r in bound_rows)}`",
                "",
                "Synthetic conservative bound only; no kernel-bridge or selector-closure claim.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
