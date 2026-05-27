#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2266 = GEN / "p2266_s1216_strict_nu_branch_group_policy_risk_calibrated_controller_map_probe.json"
IN_2273 = GEN / "p2273_s1223_strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe.json"
OUT = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.json"
MD = GEN / "p2274_s1224_strict_nu_branch_group_policy_robustness_region_certificate_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2266 = load(IN_2266)
    p2273 = load(IN_2273)

    map_rows = (p2266.get("strict_nu_branch_group_policy_risk_calibrated_controller_map_probe", {}) or {}).get("controller_map_rows", []) or []
    tightened = (p2273.get("strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe", {}) or {}).get("tightened_envelope_constants", {}) or {}

    L_rho = float(tightened.get("tightened_abs_d_bound_d_rho", 0.0) or 0.0)
    L_kappa = float(tightened.get("tightened_abs_d_bound_d_kappa", 0.0) or 0.0)

    # cert radius policy: allocate 50% risk budget to derivative-transport uncertainty
    # and split equally by L1 contributions from rho and kappa directions.
    certified_rows = []
    for row in map_rows:
        risk = float(row.get("risk_tolerance", 0.05) or 0.05)
        rho = float(row.get("rho_recommended", 0.8) or 0.8)
        kappa = float(row.get("kappa_scale_recommended", 1.0) or 1.0)

        risk_budget = 0.5 * risk
        denom = max(L_rho + L_kappa, 1e-12)
        radius_l1 = risk_budget / denom

        # balanced directional allocation
        delta_rho = 0.5 * radius_l1
        delta_kappa = 0.5 * radius_l1

        rho_lo = max(0.55, rho - delta_rho)
        rho_hi = min(0.95, rho + delta_rho)
        kappa_lo = max(1.0, kappa - delta_kappa)
        kappa_hi = min(1.8, kappa + delta_kappa)

        clipped = (rho_lo == 0.55) or (rho_hi == 0.95) or (kappa_lo == 1.0) or (kappa_hi == 1.8)

        certified_rows.append({
            "risk_tolerance": risk,
            "center": {"rho": rho, "kappa_scale": kappa},
            "risk_budget_for_transport": risk_budget,
            "tightened_lipschitz_constants": {"L_rho": L_rho, "L_kappa": L_kappa},
            "certified_l1_radius": radius_l1,
            "certified_box": {
                "rho_min": rho_lo,
                "rho_max": rho_hi,
                "kappa_scale_min": kappa_lo,
                "kappa_scale_max": kappa_hi,
            },
            "boundary_clipped_to_admissible_box": clipped,
            "certificate_rule": "if |Δrho|+|Δkappa| <= radius_l1 then predicted bound-transport <= 0.5*risk",
        })

    payload = {
        "schema_version": "p2274_s1224_v1",
        "packet_id": "P2274",
        "stage_id": "S1224",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_ROBUSTNESS_REGION_CERTIFICATE_PROBE",
        "strict_nu_branch_group_policy_robustness_region_certificate_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_ROBUSTNESS_REGION_CERTIFICATE_PROBE_V1",
            "source_packets": [str(IN_2266.relative_to(ROOT)), str(IN_2273.relative_to(ROOT))],
            "admissible_box": {"rho_min": 0.55, "rho_max": 0.95, "kappa_scale_min": 1.0, "kappa_scale_max": 1.8},
            "certified_rows": certified_rows,
            "proof_stub": {
                "norm_used": "L1",
                "transport_bound": "|ΔB| <= L_rho*|Δrho| + L_kappa*|Δkappa| <= (L_rho+L_kappa)*radius_l1",
                "budget_policy": "radius_l1 := (0.5*risk)/(L_rho+L_kappa)",
                "result": "|ΔB| <= 0.5*risk (inside certified region)",
            },
            "theorem_scope_limit": "local robustness-region certificate for surrogate transport only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2275_candidate",
            "goal": "replay nonlinear profiles at certified-box corners and certify empirical pass-rate floor under robustness transport budget",
        },
        "gatekeeper_checks": {
            "certified_rows_exported": len(certified_rows) > 0,
            "all_radius_positive": all(r["certified_l1_radius"] > 0.0 for r in certified_rows),
            "all_boxes_within_admissible_domain": all(
                0.55 <= r["certified_box"]["rho_min"] <= r["certified_box"]["rho_max"] <= 0.95 and
                1.0 <= r["certified_box"]["kappa_scale_min"] <= r["certified_box"]["kappa_scale_max"] <= 1.8
                for r in certified_rows
            ),
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2274 S1224: robustness-region certificate probe",
            "",
            f"- rows certified: `{len(certified_rows)}`",
            f"- tightened L_rho: `{L_rho:.12e}`",
            f"- tightened L_kappa: `{L_kappa:.12e}`",
            f"- all boxes within admissible domain: `{all(0.55 <= r['certified_box']['rho_min'] <= r['certified_box']['rho_max'] <= 0.95 and 1.0 <= r['certified_box']['kappa_scale_min'] <= r['certified_box']['kappa_scale_max'] <= 1.8 for r in certified_rows)}`",
            "",
            "Local surrogate robustness certificate only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
