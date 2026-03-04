#!/usr/bin/env python3
"""
QW-2094: strict-rigor defect sweep over the non-anchor closure path.

Purpose:
- verify internal consistency after QW-2093 -> QW-2085/86/87 -> QW-2084 -> QW-2069/2071,
- detect implementation/report mismatches that could create false success claims.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2094_strict_rigor_defect_sweep.json"
OUT_MD = ROOT / "RAPORT_QW2094_STRICT_RIGOR_DEFECT_SWEEP.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def bool_pass_count(flags: Dict) -> int:
    return sum(1 for v in flags.values() if bool(v))


def main() -> None:
    r2084 = load_json("report_qw2084_t1_nonanchor_strict_gate.json")
    r2085 = load_json("report_qw2085_gf_nonanchor_lifetime_gate.json")
    r2086 = load_json("report_qw2086_mz_nonanchor_ew_pole_gate.json")
    r2087 = load_json("report_qw2087_alpha_s_nonanchor_boundary_gate.json")
    r2090 = load_optional_json("report_qw2090_h0_lambda_decoupling_gate.json")
    r2091 = load_optional_json("report_qw2091_neutrino_absolute_scale_gate.json")
    r2092 = load_optional_json("report_qw2092_gnewton_si_bridge_gate.json")
    r2097 = load_optional_json("report_qw2097_ckm_cp_target_refinement_gate.json")
    r2098 = load_optional_json("report_qw2098_ew_secondary_nonanchor_closure_gate.json")
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2071 = load_json("report_qw2071_sm_gr_full_precision_closure_gate.json")

    defects: List[Dict] = []
    checks: Dict[str, bool] = {}

    # 1) Single-gate integrity checks.
    for tag, rep in [("QW-2085", r2085), ("QW-2086", r2086), ("QW-2087", r2087)]:
        flags = rep.get("flags", {})
        strict_flag = bool(flags.get("strict_nonanchor_pass", False))
        verdict_pass = str(rep.get("verdict", "")).endswith("_PASS")
        update = rep.get("update", {})
        strict_level_ok = str(update.get("strict_level", "")) == "strict_internal_gate"
        status_ok = str(update.get("status", "")).startswith("derived")
        pass_count_ok = bool_pass_count(flags) == int(rep.get("pass_count", -1))

        checks[f"{tag}_verdict_matches_flag"] = bool(verdict_pass == strict_flag)
        checks[f"{tag}_update_strict_level_ok"] = strict_level_ok
        checks[f"{tag}_update_status_ok"] = status_ok
        checks[f"{tag}_pass_count_consistent"] = pass_count_ok

    # 2) Aggregate gate integrity checks.
    flags84 = r2084.get("flags", {})
    checks["QW-2084_pass_count_consistent"] = bool_pass_count(flags84) == int(r2084.get("pass_count", -1))
    checks["QW-2084_q2093_artifacts_present"] = bool(flags84.get("q2093_kernel_input_artifacts_present", False))
    checks["QW-2084_no_anchor_or_circularity_detected"] = bool(flags84.get("no_anchor_or_circularity_detected", False))
    checks["QW-2084_alpha_s_strict_pass"] = bool(flags84.get("alpha_s_strict_nonanchor_pass", False))
    checks["QW-2084_gf_strict_pass"] = bool(flags84.get("g_f_strict_nonanchor_pass", False))
    checks["QW-2084_mz_strict_pass"] = bool(flags84.get("m_z_strict_nonanchor_pass", False))

    verdict84 = str(r2084.get("verdict", ""))
    checks["QW-2084_verdict_pass_when_all_flags"] = (
        verdict84 == "T1_NONANCHOR_STRICT_GATE_PASS"
        and checks["QW-2084_q2093_artifacts_present"]
        and checks["QW-2084_no_anchor_or_circularity_detected"]
        and checks["QW-2084_alpha_s_strict_pass"]
        and checks["QW-2084_gf_strict_pass"]
        and checks["QW-2084_mz_strict_pass"]
    )

    # 3) Cross-report consistency.
    coverage69 = r2069.get("coverage", {})
    strict_unresolved_69 = int(coverage69.get("n_strict_unresolved", -1))
    strict_unresolved_71 = len(r2071.get("strict_unresolved_parameters", []))
    checks["QW-2069_vs_QW-2071_strict_unresolved_consistent"] = strict_unresolved_69 == strict_unresolved_71

    missing_69 = int(coverage69.get("n_missing", -1))
    missing_71 = len(r2071.get("missing_parameters", []))
    checks["QW-2069_vs_QW-2071_missing_consistent"] = missing_69 == missing_71

    # 4) Ensure T1 parameters are strict in package report.
    entry_map = {e["id"]: e for e in r2069.get("entries", []) if isinstance(e, dict) and "id" in e}
    for pid in ["alpha_s_mz", "g_f", "m_z"]:
        e = entry_map.get(pid, {})
        checks[f"QW-2069_{pid}_strict_level"] = str(e.get("strict_level", "")) == "strict_internal_gate"
        checks[f"QW-2069_{pid}_status"] = str(e.get("status", "")).startswith("derived")

    # 5) Optional CKM CP target-refinement gate checks.
    if r2097 is not None:
        flags97 = r2097.get("flags", {})
        checks["QW-2097_pass_count_consistent"] = bool_pass_count(flags97) == int(
            r2097.get("pass_count", -1)
        )
        upd97 = r2097.get("update", {})
        checks["QW-2097_update_id_delta_cp_ckm"] = str(upd97.get("id", "")) == "delta_cp_ckm"
        checks["QW-2097_update_strict_level_ok"] = (
            str(upd97.get("strict_level", "")) == "strict_internal_gate"
        )
        verdict97 = str(r2097.get("verdict", ""))
        status97 = str(upd97.get("status", ""))
        if verdict97 == "CKM_CP_TARGET_REFINEMENT_GATE_PASS_STRICT":
            checks["QW-2097_verdict_status_consistent"] = status97.startswith("derived") and (
                status97 != "derived_strict_target_miss"
            )
        elif verdict97 == "CKM_CP_TARGET_REFINEMENT_GATE_TARGET_MISS":
            checks["QW-2097_verdict_status_consistent"] = status97 == "derived_strict_target_miss"
        else:
            checks["QW-2097_verdict_status_consistent"] = False

        e_ckm = entry_map.get("delta_cp_ckm", {})
        checks["QW-2069_delta_cp_ckm_matches_qw2097_status"] = (
            str(e_ckm.get("status", "")) == status97
        )
        checks["QW-2069_delta_cp_ckm_matches_qw2097_method"] = (
            str(e_ckm.get("method", "")) == str(upd97.get("method", ""))
        )

    # 6) Optional H0/Lambda decoupling gate checks (QW-2090).
    if r2090 is not None:
        flags90 = r2090.get("flags", {})
        checks["QW-2090_pass_count_consistent"] = bool_pass_count(flags90) == int(
            r2090.get("pass_count", -1)
        )
        upd90 = [u for u in r2090.get("updates", []) if isinstance(u, dict) and "id" in u]
        upd90_map = {u["id"]: u for u in upd90}
        checks["QW-2090_updates_cover_h0_lambda"] = set(upd90_map.keys()) == {
            "h0",
            "lambda_cosmological",
        }
        verdict90 = str(r2090.get("verdict", ""))
        if verdict90 == "H0_LAMBDA_DECOUPLING_GATE_PASS_STRICT":
            checks["QW-2090_verdict_status_consistent"] = all(
                str(upd90_map.get(pid, {}).get("strict_level", "")) == "strict_internal_gate"
                and str(upd90_map.get(pid, {}).get("status", "")).startswith("derived")
                and str(upd90_map.get(pid, {}).get("status", "")) != "derived_strict_target_miss"
                for pid in ["h0", "lambda_cosmological"]
            )
        elif verdict90 == "H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS":
            statuses90 = [
                str(upd90_map.get(pid, {}).get("status", ""))
                for pid in ["h0", "lambda_cosmological"]
            ]
            levels90 = [
                str(upd90_map.get(pid, {}).get("strict_level", ""))
                for pid in ["h0", "lambda_cosmological"]
            ]
            checks["QW-2090_verdict_status_consistent"] = (
                all(lvl == "strict_internal_gate" for lvl in levels90)
                and all(st.startswith("derived") for st in statuses90)
                and any(st == "derived_strict_target_miss" for st in statuses90)
            )
        elif verdict90 == "H0_LAMBDA_DECOUPLING_GATE_PENDING_NONCLOSING":
            checks["QW-2090_verdict_status_consistent"] = all(
                str(upd90_map.get(pid, {}).get("strict_level", "")) == "coupled_anchor_dependent"
                and str(upd90_map.get(pid, {}).get("status", "")) == "derived_coupled_anchor_dependent"
                for pid in ["h0", "lambda_cosmological"]
            )
        else:
            checks["QW-2090_verdict_status_consistent"] = False
        checks["QW-2069_h0_matches_qw2090_status"] = (
            str(entry_map.get("h0", {}).get("status", ""))
            == str(upd90_map.get("h0", {}).get("status", ""))
        )
        checks["QW-2069_lambda_matches_qw2090_status"] = (
            str(entry_map.get("lambda_cosmological", {}).get("status", ""))
            == str(upd90_map.get("lambda_cosmological", {}).get("status", ""))
        )

    # 7) Optional neutrino absolute-scale gate checks (QW-2091).
    if r2091 is not None:
        flags91 = r2091.get("flags", {})
        checks["QW-2091_pass_count_consistent"] = bool_pass_count(flags91) == int(
            r2091.get("pass_count", -1)
        )
        upd91 = [u for u in r2091.get("updates", []) if isinstance(u, dict) and "id" in u]
        upd91_map = {u["id"]: u for u in upd91}
        checks["QW-2091_updates_cover_neutrino_triad"] = set(upd91_map.keys()) == {
            "m_nu1",
            "m_nu2",
            "m_nu3",
        }
        verdict91 = str(r2091.get("verdict", ""))
        if verdict91 == "NEUTRINO_ABSOLUTE_SCALE_GATE_PASS_STRICT":
            checks["QW-2091_verdict_status_consistent"] = all(
                str(upd91_map.get(pid, {}).get("strict_level", "")) == "strict_internal_gate"
                and str(upd91_map.get(pid, {}).get("status", "")).startswith("derived")
                and str(upd91_map.get(pid, {}).get("status", "")) != "derived_strict_target_miss"
                for pid in ["m_nu1", "m_nu2", "m_nu3"]
            )
        elif verdict91 == "NEUTRINO_ABSOLUTE_SCALE_GATE_PENDING_NONCLOSING":
            checks["QW-2091_verdict_status_consistent"] = all(
                str(upd91_map.get(pid, {}).get("strict_level", "")) == "model_assumption_anchor"
                and str(upd91_map.get(pid, {}).get("status", "")) == "derived_model_assumption_nonclosing"
                for pid in ["m_nu1", "m_nu2", "m_nu3"]
            )
        else:
            checks["QW-2091_verdict_status_consistent"] = False
        for pid in ["m_nu1", "m_nu2", "m_nu3"]:
            checks[f"QW-2069_{pid}_matches_qw2091_status"] = (
                str(entry_map.get(pid, {}).get("status", ""))
                == str(upd91_map.get(pid, {}).get("status", ""))
            )

    # 8) Optional G_newton SI-bridge gate checks (QW-2092).
    if r2092 is not None:
        flags92 = r2092.get("flags", {})
        checks["QW-2092_pass_count_consistent"] = bool_pass_count(flags92) == int(
            r2092.get("pass_count", -1)
        )
        upd92 = r2092.get("update", {})
        checks["QW-2092_update_id_gnewton"] = str(upd92.get("id", "")) == "G_newton"
        verdict92 = str(r2092.get("verdict", ""))
        status92 = str(upd92.get("status", ""))
        level92 = str(upd92.get("strict_level", ""))
        if verdict92 == "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT":
            checks["QW-2092_verdict_status_consistent"] = (
                status92.startswith("derived")
                and status92 != "derived_strict_target_miss"
                and level92 == "strict_internal_gate"
                and bool(flags92.get("bridge_not_backsolved_from_g_si", False))
            )
        elif verdict92 == "GNEWTON_SI_BRIDGE_GATE_PENDING_NONCLOSING":
            checks["QW-2092_verdict_status_consistent"] = (
                status92 == "derived_nofit_anchor_dependent"
                and level92 == "physical_relation_anchor_dependent"
            )
        else:
            checks["QW-2092_verdict_status_consistent"] = False
        checks["QW-2069_gnewton_matches_qw2092_status"] = (
            str(entry_map.get("G_newton", {}).get("status", "")) == status92
        )

    # 9) Optional EW secondary closure gate checks (QW-2098).
    if r2098 is not None:
        flags98 = r2098.get("flags", {})
        checks["QW-2098_pass_count_consistent"] = bool_pass_count(flags98) == int(
            r2098.get("pass_count", -1)
        )
        upd98 = [u for u in r2098.get("updates", []) if isinstance(u, dict) and "id" in u]
        upd98_map = {u["id"]: u for u in upd98}
        exp_98 = {"v_higgs", "m_w", "sin2_theta_w_mz", "alpha_em_inv_mz"}
        checks["QW-2098_updates_cover_ew_secondary_set"] = set(upd98_map.keys()) == exp_98

        verdict98 = str(r2098.get("verdict", ""))
        statuses98 = [str(upd98_map.get(pid, {}).get("status", "")) for pid in exp_98]
        levels98 = [str(upd98_map.get(pid, {}).get("strict_level", "")) for pid in exp_98]
        if verdict98 == "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_PASS_STRICT":
            checks["QW-2098_verdict_status_consistent"] = all(
                s.startswith("derived") and s != "derived_strict_target_miss" for s in statuses98
            ) and all(l == "strict_internal_gate" for l in levels98)
        elif verdict98 == "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_TARGET_MISS":
            checks["QW-2098_verdict_status_consistent"] = (
                all(l == "strict_internal_gate" for l in levels98)
                and all(s.startswith("derived") for s in statuses98)
                and any(s == "derived_strict_target_miss" for s in statuses98)
            )
        elif verdict98 == "EW_SECONDARY_NONANCHOR_CLOSURE_GATE_FAIL_CHAIN_NOT_READY":
            checks["QW-2098_verdict_status_consistent"] = (
                all(l == "physical_relation_anchor_dependent" for l in levels98)
                and all(s == "derived_nofit_anchor_dependent" for s in statuses98)
            )
        else:
            checks["QW-2098_verdict_status_consistent"] = False

        for pid in exp_98:
            checks[f"QW-2069_{pid}_matches_qw2098_status"] = (
                str(entry_map.get(pid, {}).get("status", ""))
                == str(upd98_map.get(pid, {}).get("status", ""))
            )
            checks[f"QW-2069_{pid}_matches_qw2098_method"] = (
                str(entry_map.get(pid, {}).get("method", ""))
                == str(upd98_map.get(pid, {}).get("method", ""))
            )

    # Build defect list.
    for name, ok in checks.items():
        if ok:
            continue
        defects.append(
            {
                "id": name,
                "severity": "high",
                "message": "Consistency check failed; potential false-success risk.",
            }
        )

    verdict = (
        "STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS"
        if not defects
        else "STRICT_RIGOR_DEFECT_SWEEP_FAIL_CRITICAL_DEFECTS_FOUND"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2084": "report_qw2084_t1_nonanchor_strict_gate.json",
            "qw2085": "report_qw2085_gf_nonanchor_lifetime_gate.json",
            "qw2086": "report_qw2086_mz_nonanchor_ew_pole_gate.json",
            "qw2087": "report_qw2087_alpha_s_nonanchor_boundary_gate.json",
            "qw2090": (
                "report_qw2090_h0_lambda_decoupling_gate.json"
                if r2090 is not None
                else None
            ),
            "qw2091": (
                "report_qw2091_neutrino_absolute_scale_gate.json"
                if r2091 is not None
                else None
            ),
            "qw2092": (
                "report_qw2092_gnewton_si_bridge_gate.json"
                if r2092 is not None
                else None
            ),
            "qw2097": (
                "report_qw2097_ckm_cp_target_refinement_gate.json"
                if r2097 is not None
                else None
            ),
            "qw2098": (
                "report_qw2098_ew_secondary_nonanchor_closure_gate.json"
                if r2098 is not None
                else None
            ),
            "qw2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "qw2071": "report_qw2071_sm_gr_full_precision_closure_gate.json",
        },
        "checks": checks,
        "n_checks": len(checks),
        "n_failed_checks": len(defects),
        "defects": defects,
        "verdict": verdict,
        "required_next_step": (
            "NO_CRITICAL_DEFECTS_IN_THIS_SWEEP"
            if verdict == "STRICT_RIGOR_DEFECT_SWEEP_PASS_NO_CRITICAL_DEFECTS"
            else "FIX_FAILED_CHECKS_AND_RERUN_SWEEP"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2094: STRICT RIGOR DEFECT SWEEP",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Checks: `{out['n_checks']}`",
        f"- Failed checks: `{out['n_failed_checks']}`",
        "",
        "## Failed Checks",
    ]
    if defects:
        for d in defects:
            lines.append(f"- `{d['id']}` ({d['severity']}): {d['message']}")
    else:
        lines.append("- none")
    lines += [
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2094] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2094] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2094] verdict={verdict} checks={out['n_checks']} failed={out['n_failed_checks']}")


if __name__ == "__main__":
    main()
