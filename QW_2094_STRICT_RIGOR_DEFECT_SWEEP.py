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
    r2099 = load_optional_json("report_qw2099_hz_external_decoupling_autocollector.json")
    r2101 = load_optional_json("report_qw2101_gnewton_bridge_external_autocollector.json")
    r2097 = load_optional_json("report_qw2097_ckm_cp_target_refinement_gate.json")
    r2098 = load_optional_json("report_qw2098_ew_secondary_nonanchor_closure_gate.json")
    r2102 = load_optional_json("report_qw2102_hz_decoupling_identifiability_gate.json")
    r2103 = load_optional_json("report_qw2103_gnewton_dimensionless_provenance_gate.json")
    r2104 = load_optional_json("report_qw2104_t3t4_strict_preflight_gate.json")
    r2105 = load_optional_json("report_qw2105_t3t4_strict_input_gap_report.json")
    r2106 = load_optional_json("report_qw2106_strict_external_input_intake_gate.json")
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

    # 10) Optional H(z) identifiability pre-gate checks (QW-2102).
    if r2102 is not None:
        flags102 = r2102.get("flags", {})
        met102 = r2102.get("metrics", {})
        checks["QW-2102_pass_count_consistent"] = bool_pass_count(flags102) == int(
            r2102.get("pass_count", -1)
        )

        verdict102 = str(r2102.get("verdict", ""))
        checks["QW-2102_verdict_known"] = verdict102 in {
            "HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY",
            "HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING",
        }
        if verdict102 == "HZ_DECOUPLING_IDENTIFIABILITY_GATE_WEAK_LEVERARM_PENDING":
            checks["QW-2102_weak_blocks_qw2090_strict_pass"] = (
                r2090 is None
                or str(r2090.get("verdict", "")) != "H0_LAMBDA_DECOUPLING_GATE_PASS_STRICT"
            )
        elif verdict102 == "HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY":
            checks["QW-2102_weak_blocks_qw2090_strict_pass"] = True
        else:
            checks["QW-2102_weak_blocks_qw2090_strict_pass"] = False

        n_nodes = float(met102.get("n_nodes", -1.0))
        z_span = met102.get("z_span")
        e_span = met102.get("e_span")
        cond = met102.get("design_condition_number")
        checks["QW-2102_nodes_flag_matches_metric"] = (
            bool(flags102.get("n_nodes_ge_5", False)) == bool(n_nodes >= 5.0)
        )
        checks["QW-2102_zspan_flag_matches_metric"] = (
            bool(flags102.get("z_span_ge_0p8", False))
            == bool(isinstance(z_span, (int, float)) and float(z_span) >= 0.8)
        )
        checks["QW-2102_espan_flag_matches_metric"] = (
            bool(flags102.get("e_span_ge_1p0", False))
            == bool(isinstance(e_span, (int, float)) and float(e_span) >= 1.0)
        )
        checks["QW-2102_cond_flag_matches_metric"] = (
            bool(flags102.get("design_condition_lt_8", False))
            == bool(isinstance(cond, (int, float)) and float(cond) < 8.0)
        )
        ctx102 = r2102.get("qw2090_context", {})
        checks["QW-2102_qw2090_context_verdict_consistent"] = (
            r2090 is None
            or str(ctx102.get("qw2090_verdict", "")) == str(r2090.get("verdict", ""))
        )

    # 11) Optional G_newton provenance pre-gate checks (QW-2103).
    if r2103 is not None:
        flags103 = r2103.get("flags", {})
        checks["QW-2103_pass_count_consistent"] = bool_pass_count(flags103) == int(
            r2103.get("pass_count", -1)
        )

        verdict103 = str(r2103.get("verdict", ""))
        checks["QW-2103_verdict_known"] = verdict103 in {
            "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY",
            "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING",
        }
        if verdict103 == "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING":
            checks["QW-2103_pending_blocks_qw2092_strict_pass"] = (
                r2092 is None
                or str(r2092.get("verdict", "")) != "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT"
            )
        elif verdict103 == "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY":
            checks["QW-2103_pending_blocks_qw2092_strict_pass"] = True
        else:
            checks["QW-2103_pending_blocks_qw2092_strict_pass"] = False

        summary103 = r2103.get("input_summary", {})
        origin103 = str(summary103.get("bridge_observable_origin", ""))
        checks["QW-2103_origin_matches_flag"] = (
            (origin103 == "external_dimensionless_observable")
            == bool(flags103.get("bridge_origin_external_dimensionless", False))
        )
        checks["QW-2103_anchor_free_matches_flag"] = (
            bool(summary103.get("provenance_anchor_free", False))
            == bool(flags103.get("provenance_anchor_free", False))
        )
        gsi_opt = summary103.get("g_si_input_optional", None)
        checks["QW-2103_gsi_optional_matches_flag"] = (
            (gsi_opt is None) == bool(flags103.get("g_si_not_primary_input", False))
        )
        if r2092 is not None:
            flags92 = r2092.get("flags", {})
            checks["QW-2103_matches_qw2092_backsolve_flag"] = (
                (origin103 == "backsolved_from_g_si")
                == (not bool(flags92.get("bridge_not_backsolved_from_g_si", False)))
            )
        else:
            checks["QW-2103_matches_qw2092_backsolve_flag"] = True

    # 12) Optional merged T3/T4 preflight meta-gate checks (QW-2104).
    if r2104 is not None:
        flags104 = r2104.get("flags", {})
        checks["QW-2104_pass_count_consistent"] = bool_pass_count(flags104) == int(
            r2104.get("pass_count", -1)
        )
        verdict104 = str(r2104.get("verdict", ""))
        checks["QW-2104_verdict_known"] = verdict104 in {
            "T3T4_STRICT_PREFLIGHT_GATE_PASS",
            "T3T4_STRICT_PREFLIGHT_GATE_PENDING",
            "T3T4_STRICT_PREFLIGHT_GATE_FAIL_LOGIC_DEFECT",
        }
        defects104 = r2104.get("defects", [])
        checks["QW-2104_defects_list_type"] = isinstance(defects104, list)

        if r2099 is not None:
            checks["QW-2104_hz_input_ready_matches_qw2099"] = (
                bool(flags104.get("hz_input_strict_ready", False))
                == bool(r2099.get("strict_ready", False))
            )
        else:
            checks["QW-2104_hz_input_ready_matches_qw2099"] = True

        if r2102 is not None:
            checks["QW-2104_hz_identifiability_matches_qw2102"] = (
                bool(flags104.get("hz_identifiability_gate_pass", False))
                == (
                    str(r2102.get("verdict", ""))
                    == "HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY"
                )
            )
        else:
            checks["QW-2104_hz_identifiability_matches_qw2102"] = True

        if r2090 is not None:
            checks["QW-2104_hz_decoupling_matches_qw2090"] = (
                bool(flags104.get("hz_decoupling_gate_strict_pass", False))
                == (str(r2090.get("verdict", "")) == "H0_LAMBDA_DECOUPLING_GATE_PASS_STRICT")
            )
        else:
            checks["QW-2104_hz_decoupling_matches_qw2090"] = True

        if r2101 is not None:
            checks["QW-2104_g_bridge_ready_matches_qw2101"] = (
                bool(flags104.get("g_bridge_input_strict_ready", False))
                == bool(r2101.get("bridge", {}).get("strict_provenance_ready", False))
            )
        else:
            checks["QW-2104_g_bridge_ready_matches_qw2101"] = True

        if r2103 is not None:
            checks["QW-2104_g_provenance_matches_qw2103"] = (
                bool(flags104.get("g_provenance_gate_pass", False))
                == (
                    str(r2103.get("verdict", ""))
                    == "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY"
                )
            )
        else:
            checks["QW-2104_g_provenance_matches_qw2103"] = True

        if r2092 is not None:
            checks["QW-2104_g_si_bridge_matches_qw2092"] = (
                bool(flags104.get("g_si_bridge_gate_strict_pass", False))
                == (str(r2092.get("verdict", "")) == "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT")
            )
        else:
            checks["QW-2104_g_si_bridge_matches_qw2092"] = True

        all104 = bool(all(bool(v) for v in flags104.values()))
        if verdict104 == "T3T4_STRICT_PREFLIGHT_GATE_PASS":
            checks["QW-2104_verdict_consistent_with_flags"] = all104 and len(defects104) == 0
        elif verdict104 == "T3T4_STRICT_PREFLIGHT_GATE_PENDING":
            checks["QW-2104_verdict_consistent_with_flags"] = (not all104) and len(defects104) == 0
        elif verdict104 == "T3T4_STRICT_PREFLIGHT_GATE_FAIL_LOGIC_DEFECT":
            checks["QW-2104_verdict_consistent_with_flags"] = len(defects104) > 0
        else:
            checks["QW-2104_verdict_consistent_with_flags"] = False

    # 13) Optional strict external raw-input intake checks (QW-2106).
    if r2106 is not None:
        hz106 = r2106.get("hz_flags", {})
        g106 = r2106.get("g_flags", {})
        all106 = {**hz106, **g106}
        checks["QW-2106_pass_count_consistent"] = bool_pass_count(all106) == int(
            r2106.get("pass_count", -1)
        )
        verdict106 = str(r2106.get("verdict", ""))
        checks["QW-2106_verdict_known"] = verdict106 in {
            "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS",
            "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING",
        }
        all_flags_true106 = all(bool(v) for v in all106.values()) if all106 else False
        if verdict106 == "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS":
            checks["QW-2106_verdict_consistent_with_flags"] = all_flags_true106
        elif verdict106 == "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING":
            checks["QW-2106_verdict_consistent_with_flags"] = not all_flags_true106
        else:
            checks["QW-2106_verdict_consistent_with_flags"] = False

        if r2102 is not None:
            f102 = r2102.get("flags", {})
            checks["QW-2106_hz_nodes_flag_matches_qw2102"] = (
                bool(hz106.get("hz_n_nodes_ge_5", False)) == bool(f102.get("n_nodes_ge_5", False))
            )
            checks["QW-2106_hz_zspan_flag_matches_qw2102"] = (
                bool(hz106.get("hz_z_span_ge_0p8", False)) == bool(f102.get("z_span_ge_0p8", False))
            )
            checks["QW-2106_hz_espan_flag_matches_qw2102"] = (
                bool(hz106.get("hz_e_span_ge_1p0", False)) == bool(f102.get("e_span_ge_1p0", False))
            )
            checks["QW-2106_hz_cond_flag_matches_qw2102"] = (
                bool(hz106.get("hz_design_condition_lt_8", False))
                == bool(f102.get("design_condition_lt_8", False))
            )
        else:
            checks["QW-2106_hz_nodes_flag_matches_qw2102"] = True
            checks["QW-2106_hz_zspan_flag_matches_qw2102"] = True
            checks["QW-2106_hz_espan_flag_matches_qw2102"] = True
            checks["QW-2106_hz_cond_flag_matches_qw2102"] = True

        if r2103 is not None:
            f103 = r2103.get("flags", {})
            checks["QW-2106_g_origin_flag_matches_qw2103"] = (
                bool(g106.get("g_origin_external_dimensionless", False))
                == bool(f103.get("bridge_origin_external_dimensionless", False))
            )
            checks["QW-2106_g_anchor_free_flag_matches_qw2103"] = (
                bool(g106.get("g_provenance_anchor_free", False))
                == bool(f103.get("provenance_anchor_free", False))
            )
            checks["QW-2106_g_not_seeded_flag_matches_qw2103"] = (
                bool(g106.get("g_not_seeded_from_registry", False))
                == bool(f103.get("not_seeded_from_registry", False))
            )
            checks["QW-2106_g_si_primary_flag_matches_qw2103"] = (
                bool(g106.get("g_si_not_primary", False)) == bool(f103.get("g_si_not_primary_input", False))
            )
        else:
            checks["QW-2106_g_origin_flag_matches_qw2103"] = True
            checks["QW-2106_g_anchor_free_flag_matches_qw2103"] = True
            checks["QW-2106_g_not_seeded_flag_matches_qw2103"] = True
            checks["QW-2106_g_si_primary_flag_matches_qw2103"] = True

        if r2104 is not None:
            verdict104 = str(r2104.get("verdict", ""))
            checks["QW-2106_pending_blocks_qw2104_pass"] = not (
                verdict106 == "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PENDING"
                and verdict104 == "T3T4_STRICT_PREFLIGHT_GATE_PASS"
            )
        else:
            checks["QW-2106_pending_blocks_qw2104_pass"] = True

    # 14) Optional T3/T4 gap-report consistency checks (QW-2105).
    if r2105 is not None:
        verdict105 = str(r2105.get("verdict", ""))
        checks["QW-2105_verdict_known"] = verdict105 in {
            "T3T4_STRICT_INPUT_GAPS_PRESENT",
            "T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN",
        }
        hz_path = r2105.get("hz_path", {})
        g_path = r2105.get("gnewton_path", {})
        hz_gaps = hz_path.get("gaps", [])
        g_gaps = g_path.get("gaps", [])
        hz_ready = bool(hz_path.get("strict_ready", False))
        g_ready = bool(g_path.get("strict_ready", False))
        checks["QW-2105_hz_ready_matches_gaps"] = isinstance(hz_gaps, list) and (hz_ready == (len(hz_gaps) == 0))
        checks["QW-2105_g_ready_matches_gaps"] = isinstance(g_gaps, list) and (g_ready == (len(g_gaps) == 0))

        if verdict105 == "T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN":
            checks["QW-2105_verdict_consistent_with_ready"] = hz_ready and g_ready
        elif verdict105 == "T3T4_STRICT_INPUT_GAPS_PRESENT":
            checks["QW-2105_verdict_consistent_with_ready"] = not (hz_ready and g_ready)
        else:
            checks["QW-2105_verdict_consistent_with_ready"] = False

        if r2104 is not None:
            verdict104 = str(r2104.get("verdict", ""))
            checks["QW-2105_gaps_present_blocks_qw2104_pass"] = not (
                verdict105 == "T3T4_STRICT_INPUT_GAPS_PRESENT"
                and verdict104 == "T3T4_STRICT_PREFLIGHT_GATE_PASS"
            )
        else:
            checks["QW-2105_gaps_present_blocks_qw2104_pass"] = True

        if r2106 is not None:
            intake106 = r2106.get("verdict", "")
            checks["QW-2105_closed_implies_qw2106_pass"] = not (
                verdict105 == "T3T4_STRICT_INPUT_GAPS_CLOSED_READY_FOR_STRICT_RERUN"
                and intake106 != "STRICT_EXTERNAL_INPUT_INTAKE_GATE_PASS"
            )
        else:
            checks["QW-2105_closed_implies_qw2106_pass"] = True

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
            "qw2099": (
                "report_qw2099_hz_external_decoupling_autocollector.json"
                if r2099 is not None
                else None
            ),
            "qw2101": (
                "report_qw2101_gnewton_bridge_external_autocollector.json"
                if r2101 is not None
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
            "qw2102": (
                "report_qw2102_hz_decoupling_identifiability_gate.json"
                if r2102 is not None
                else None
            ),
            "qw2103": (
                "report_qw2103_gnewton_dimensionless_provenance_gate.json"
                if r2103 is not None
                else None
            ),
            "qw2104": (
                "report_qw2104_t3t4_strict_preflight_gate.json"
                if r2104 is not None
                else None
            ),
            "qw2105": (
                "report_qw2105_t3t4_strict_input_gap_report.json"
                if r2105 is not None
                else None
            ),
            "qw2106": (
                "report_qw2106_strict_external_input_intake_gate.json"
                if r2106 is not None
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
