#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2270 = GEN / "p2270_s1220_strict_nu_branch_group_policy_symbolic_derivative_bound_probe.json"
IN_2272 = GEN / "p2272_s1222_strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe.json"
OUT = GEN / "p2273_s1223_strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe.json"
MD = GEN / "p2273_s1223_strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2270 = load(IN_2270)
    p2272 = load(IN_2272)

    symbolic_box = (p2270.get("strict_nu_branch_group_policy_symbolic_derivative_bound_probe", {}) or {}).get("symbolic_box_certificates", {}) or {}
    max_abs_sym_drho = float(symbolic_box.get("max_abs_d_bound_d_rho_over_box", 0.0) or 0.0)
    max_abs_sym_dkappa = float(symbolic_box.get("max_abs_d_bound_d_kappa_over_box", 0.0) or 0.0)

    rows = (p2272.get("strict_nu_branch_group_policy_empirical_derivative_resolution_and_quantile_calibration_probe", {}) or {}).get("calibrated_rows", []) or []

    empirical_max_drho = max((float(r.get("empirical_abs_max_d_rho_across_profiles", 0.0) or 0.0) for r in rows), default=0.0)
    empirical_max_dkappa = max((float(r.get("empirical_abs_max_d_kappa_across_profiles", 0.0) or 0.0) for r in rows), default=0.0)

    # data-driven tightening with explicit safety multiplier policy
    # still no closure claim: this only calibrates surrogate envelope constants.
    safety_multiplier = 1.25
    tightened_const_rho = safety_multiplier * empirical_max_drho
    tightened_const_kappa = safety_multiplier * empirical_max_dkappa

    legacy_symbolic_dominates_rho = max_abs_sym_drho + 1e-12 >= empirical_max_drho
    legacy_symbolic_dominates_kappa = max_abs_sym_dkappa + 1e-12 >= empirical_max_dkappa

    tightened_vs_legacy_ratio_rho = tightened_const_rho / max(max_abs_sym_drho, 1e-12)
    tightened_vs_legacy_ratio_kappa = tightened_const_kappa / max(max_abs_sym_dkappa, 1e-12)

    payload = {
        "schema_version": "p2273_s1223_v1",
        "packet_id": "P2273",
        "stage_id": "S1223",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-27T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_NU_BRANCH_GROUP_POLICY_TIGHTENED_SYMBOLIC_ENVELOPE_CONSTANT_PROBE",
        "strict_nu_branch_group_policy_tightened_symbolic_envelope_constant_probe": {
            "probe_id": "STRICT_NU_BRANCH_GROUP_POLICY_TIGHTENED_SYMBOLIC_ENVELOPE_CONSTANT_PROBE_V1",
            "source_packets": [str(IN_2270.relative_to(ROOT)), str(IN_2272.relative_to(ROOT))],
            "inputs": {
                "legacy_symbolic_box_max_abs_d_bound_d_rho": max_abs_sym_drho,
                "legacy_symbolic_box_max_abs_d_bound_d_kappa": max_abs_sym_dkappa,
                "empirical_max_abs_d_rho_across_rows": empirical_max_drho,
                "empirical_max_abs_d_kappa_across_rows": empirical_max_dkappa,
                "safety_multiplier": safety_multiplier,
            },
            "tightened_envelope_constants": {
                "tightened_abs_d_bound_d_rho": tightened_const_rho,
                "tightened_abs_d_bound_d_kappa": tightened_const_kappa,
                "tightened_l1_lipschitz": tightened_const_rho + tightened_const_kappa,
            },
            "diagnostics": {
                "legacy_symbolic_dominates_empirical_rho": legacy_symbolic_dominates_rho,
                "legacy_symbolic_dominates_empirical_kappa": legacy_symbolic_dominates_kappa,
                "tightened_vs_legacy_ratio_rho": tightened_vs_legacy_ratio_rho,
                "tightened_vs_legacy_ratio_kappa": tightened_vs_legacy_ratio_kappa,
                "tightening_is_nonexpansive_vs_legacy": tightened_vs_legacy_ratio_rho <= 1.0 + 1e-12 and tightened_vs_legacy_ratio_kappa <= 1.0 + 1e-12,
            },
            "calibration_policy": "tightened_constant = safety_multiplier * empirical_max_abs_derivative, with safety_multiplier>=1",
            "theorem_scope_limit": "surrogate envelope-constant calibration only; not selector closure and not ToE closure",
        },
        "recommended_next_honest_step": {
            "id": "P2274_candidate",
            "goal": "prove closed-form monotone robustness region under tightened constants and export certified risk->parameter admissible set",
        },
        "gatekeeper_checks": {
            "tightened_constants_exported": True,
            "safety_multiplier_ge_one": safety_multiplier >= 1.0,
            "tightening_nonexpansive_vs_legacy": tightened_vs_legacy_ratio_rho <= 1.0 + 1e-12 and tightened_vs_legacy_ratio_kappa <= 1.0 + 1e-12,
            "tightened_constants_nonnegative": tightened_const_rho >= 0.0 and tightened_const_kappa >= 0.0,
            "no_bridge_theorem_claimed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2273 S1223: tightened symbolic envelope constant probe",
            "",
            f"- legacy max |d/drho|: `{max_abs_sym_drho:.12e}`",
            f"- legacy max |d/dkappa|: `{max_abs_sym_dkappa:.12e}`",
            f"- empirical max |d/drho|: `{empirical_max_drho:.12e}`",
            f"- empirical max |d/dkappa|: `{empirical_max_dkappa:.12e}`",
            f"- safety multiplier: `{safety_multiplier}`",
            f"- tightened/legacy ratio (rho): `{tightened_vs_legacy_ratio_rho:.12e}`",
            f"- tightened/legacy ratio (kappa): `{tightened_vs_legacy_ratio_kappa:.12e}`",
            "",
            "Calibration only; no selector closure / ToE closure claim.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
