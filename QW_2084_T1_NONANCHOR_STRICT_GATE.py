#!/usr/bin/env python3
"""
QW-2084: T1 non-anchor strict gate (alpha_s_mz, G_F, M_Z).

Purpose:
- test whether T1 parameters can be promoted to strict non-anchor closure
  with currently available non-anchor artifacts,
- explicitly detect anchor leakage / circular dependencies,
- keep deterministic no-retune/no-scan behavior and produce auditable updates.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2084_t1_nonanchor_strict_gate.json"
OUT_MD = ROOT / "RAPORT_QW2084_T1_NONANCHOR_STRICT_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def get_registry_item(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-15) * 100.0


def _safe_float(x) -> float | None:
    if isinstance(x, (int, float)):
        return float(x)
    return None


def _extract_from_single_gate(
    pid: str,
    gate_report: Dict | None,
    ref: float,
    tol: float,
    fallback_candidate: float,
    fallback_method: str,
    fallback_reason: str,
) -> tuple[Dict, Dict]:
    if gate_report is not None:
        upd = gate_report.get("update")
        if isinstance(upd, dict) and upd.get("id") == pid:
            checks = gate_report.get("checks", {}) if isinstance(gate_report.get("checks"), dict) else {}
            flags = gate_report.get("flags", {}) if isinstance(gate_report.get("flags"), dict) else {}

            candidate = _safe_float(upd.get("predicted_value"))
            if candidate is None:
                candidate = _safe_float(checks.get(f"{pid}_candidate"))
            if candidate is None:
                candidate = ref

            rel = _safe_float(checks.get("rel_err_pct"))
            if rel is None:
                rel = rel_err_pct(candidate, ref)

            strict_pass = bool(flags.get("strict_nonanchor_pass", False))
            anchor_or_circular = bool(flags.get("anchor_leak_detected", False) or flags.get("circularity_detected", False))

            check = {
                "candidate_value": candidate,
                "registry_reference": ref,
                "rel_err_pct": rel,
                "tolerance_rel_pct": tol,
                "source_gate_verdict": gate_report.get("verdict"),
                "source_gate_pass_count": gate_report.get("pass_count"),
                "source_gate_total_flags": gate_report.get("total_flags"),
                "anchor_or_circularity_detected": anchor_or_circular,
                "strict_nonanchor_pass": strict_pass,
                "reason_if_fail": (
                    ""
                    if strict_pass
                    else "Upstream single-parameter non-anchor gate did not pass strict criteria."
                ),
            }

            update = {
                "id": pid,
                "predicted_value": candidate,
                "method": upd.get("method", f"qw2084_aggregate_from_{pid}_single_gate"),
                "status": upd.get("status", "derived_nofit_anchor_dependent"),
                "strict_level": upd.get("strict_level", "physical_relation_anchor_dependent"),
                "notes": upd.get(
                    "notes",
                    "Propagated from dedicated non-anchor gate.",
                ),
            }
            return check, update

    rel = rel_err_pct(fallback_candidate, ref)
    check = {
        "candidate_value": fallback_candidate,
        "registry_reference": ref,
        "rel_err_pct": rel,
        "tolerance_rel_pct": tol,
        "source_gate_verdict": None,
        "source_gate_pass_count": None,
        "source_gate_total_flags": None,
        "anchor_or_circularity_detected": True,
        "strict_nonanchor_pass": False,
        "reason_if_fail": fallback_reason,
    }
    update = {
        "id": pid,
        "predicted_value": fallback_candidate,
        "method": fallback_method,
        "status": "derived_nofit_anchor_dependent",
        "strict_level": "physical_relation_anchor_dependent",
        "notes": fallback_reason,
    }
    return check, update


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2070 = load_optional_json("report_qw2070_full_radiative_program_baseline.json")
    r2072 = load_optional_json("report_qw2072_ew_yukawa_flavor_radiative_baselines.json")
    r2085 = load_optional_json("report_qw2085_gf_nonanchor_lifetime_gate.json")
    r2086 = load_optional_json("report_qw2086_mz_nonanchor_ew_pole_gate.json")
    r2087 = load_optional_json("report_qw2087_alpha_s_nonanchor_boundary_gate.json")
    r2093 = load_optional_json("report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json")

    groups = reg["groups"]
    alpha_s_ref = float(get_registry_item(groups, "alpha_s_mz")["value"])
    alpha_s_tol = float(get_registry_item(groups, "alpha_s_mz")["tolerance_rel_pct"])
    gf_ref = float(get_registry_item(groups, "g_f")["value"])
    gf_tol = float(get_registry_item(groups, "g_f")["tolerance_rel_pct"])
    mz_ref = float(get_registry_item(groups, "m_z")["value"])
    mz_tol = float(get_registry_item(groups, "m_z")["tolerance_rel_pct"])

    # Fallback candidates use anchored baselines only when dedicated gates are absent.
    alpha_s_fallback = alpha_s_ref
    if r2070 is not None:
        qcd = r2070.get("qcd_one_loop_baseline", {})
        samples = qcd.get("samples", [])
        for row in samples:
            if abs(float(row.get("mu_gev", -1.0)) - mz_ref) < 1e-12:
                alpha_s_fallback = float(row["alpha_s"])
                break

    gf_fallback = gf_ref
    mz_fallback = mz_ref
    if r2072 is not None:
        ew = r2072.get("electroweak_baseline", {})
        v_from_gf = _safe_float(ew.get("v_from_gf"))
        if isinstance(v_from_gf, float) and v_from_gf > 0.0:
            gf_fallback = float(1.0 / (math.sqrt(2.0) * v_from_gf * v_from_gf))
        mw_deltar0 = _safe_float(ew.get("m_w_delta_r0_from_alpha_gf_mz"))
        sin2_in = _safe_float(ew.get("sin2_theta_w_input"))
        if isinstance(mw_deltar0, float) and isinstance(sin2_in, float) and 0.0 < sin2_in < 1.0:
            mz_fallback = float(mw_deltar0 / math.sqrt(max(1.0 - sin2_in, 1e-15)))

    alpha_s_check, alpha_s_update = _extract_from_single_gate(
        pid="alpha_s_mz",
        gate_report=r2087,
        ref=alpha_s_ref,
        tol=alpha_s_tol,
        fallback_candidate=alpha_s_fallback,
        fallback_method="qw2084_fallback_from_qw2070_qcd_baseline",
        fallback_reason="Dedicated non-anchor alpha_s gate report is unavailable or invalid.",
    )
    gf_check, gf_update = _extract_from_single_gate(
        pid="g_f",
        gate_report=r2085,
        ref=gf_ref,
        tol=gf_tol,
        fallback_candidate=gf_fallback,
        fallback_method="qw2084_fallback_gf_from_qw2072_v_identity",
        fallback_reason="Dedicated non-anchor G_F gate report is unavailable or invalid.",
    )
    mz_check, mz_update = _extract_from_single_gate(
        pid="m_z",
        gate_report=r2086,
        ref=mz_ref,
        tol=mz_tol,
        fallback_candidate=mz_fallback,
        fallback_method="qw2084_fallback_mz_from_qw2072_ew_baseline",
        fallback_reason="Dedicated non-anchor M_Z gate report is unavailable or invalid.",
    )

    alpha_s_strict_nonanchor_pass = bool(alpha_s_check["strict_nonanchor_pass"])
    gf_strict_nonanchor_pass = bool(gf_check["strict_nonanchor_pass"])
    mz_strict_nonanchor_pass = bool(mz_check["strict_nonanchor_pass"])

    checks = {
        "alpha_s_mz": alpha_s_check,
        "g_f": gf_check,
        "m_z": mz_check,
    }

    updates: List[Dict] = [alpha_s_update, gf_update, mz_update]

    q2093_artifacts_present = bool(
        r2093 is not None
        and str(r2093.get("verdict", "")) == "KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN"
        and (ROOT / "t1_nonanchor_observables_input_qw2085_2086.json").exists()
        and (ROOT / "t1_nonanchor_alpha_s_input_qw2087.json").exists()
    )

    no_anchor_or_circularity_detected = bool(
        not (
            alpha_s_check.get("anchor_or_circularity_detected", False)
            or gf_check.get("anchor_or_circularity_detected", False)
            or mz_check.get("anchor_or_circularity_detected", False)
        )
    )

    flags = {
        "deterministic_no_retune_no_scan": True,
        "q2093_kernel_input_artifacts_present": q2093_artifacts_present,
        "alpha_s_strict_nonanchor_pass": alpha_s_strict_nonanchor_pass,
        "g_f_strict_nonanchor_pass": gf_strict_nonanchor_pass,
        "m_z_strict_nonanchor_pass": mz_strict_nonanchor_pass,
        "no_anchor_or_circularity_detected": no_anchor_or_circularity_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    if (
        alpha_s_strict_nonanchor_pass
        and gf_strict_nonanchor_pass
        and mz_strict_nonanchor_pass
        and q2093_artifacts_present
        and no_anchor_or_circularity_detected
    ):
        verdict = "T1_NONANCHOR_STRICT_GATE_PASS"
    else:
        verdict = "T1_NONANCHOR_STRICT_GATE_FAIL_ANCHOR_LEAKAGE_OR_CIRCULARITY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "qcd_baseline_fallback_optional": "report_qw2070_full_radiative_program_baseline.json",
            "ew_baseline_fallback_optional": "report_qw2072_ew_yukawa_flavor_radiative_baselines.json",
            "gf_nonanchor_lifetime_gate_optional": "report_qw2085_gf_nonanchor_lifetime_gate.json",
            "mz_nonanchor_ew_pole_gate_optional": "report_qw2086_mz_nonanchor_ew_pole_gate.json",
            "alpha_s_nonanchor_boundary_gate_optional": "report_qw2087_alpha_s_nonanchor_boundary_gate.json",
            "kernel_nonanchor_input_executor_optional": "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json",
        },
        "checks": checks,
        "updates": updates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_T1_NONANCHOR_UPDATES_TO_QW2069_AND_RERUN_QW2071"
            if verdict == "T1_NONANCHOR_STRICT_GATE_PASS"
            else "PROVIDE_VALID_UPSTREAM_NONANCHOR_GATE_REPORTS_QW2085_QW2086_QW2087"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2084: T1 NONANCHOR STRICT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## T1 Checks",
        (
            f"- alpha_s_mz: strict_nonanchor_pass={alpha_s_strict_nonanchor_pass}, "
            f"rel_err_pct={float(alpha_s_check['rel_err_pct']):.6f}, "
            f"anchor_or_circular={bool(alpha_s_check['anchor_or_circularity_detected'])}"
        ),
        (
            f"- g_f: strict_nonanchor_pass={gf_strict_nonanchor_pass}, "
            f"rel_err_pct={float(gf_check['rel_err_pct']):.6f}, "
            f"anchor_or_circular={bool(gf_check['anchor_or_circularity_detected'])}"
        ),
        (
            f"- m_z: strict_nonanchor_pass={mz_strict_nonanchor_pass}, "
            f"rel_err_pct={float(mz_check['rel_err_pct']):.6f}, "
            f"anchor_or_circular={bool(mz_check['anchor_or_circularity_detected'])}"
        ),
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2084] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2084] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2084] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"alpha_s_pass={alpha_s_strict_nonanchor_pass} g_f_pass={gf_strict_nonanchor_pass} m_z_pass={mz_strict_nonanchor_pass}"
    )


if __name__ == "__main__":
    main()
