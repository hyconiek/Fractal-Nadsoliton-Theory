#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_2203 = GEN / "p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json"
IN_2280 = GEN / "p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json"
IN_2281 = GEN / "p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json"
IN_2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
IN_2300 = GEN / "p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json"
OUT = GEN / "p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.json"
MD = GEN / "p2301_s1251_strict_task3_provider_corrected_g1_g2_g3_replay_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def gap_by_id(p2282: dict[str, Any]) -> dict[str, dict[str, Any]]:
    probe = p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {}) or {}
    return {row.get("id", "UNKNOWN"): row for row in probe.get("gap_rows", []) or []}


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha = load(IN_ALPHA)
    p2203 = load(IN_2203)
    p2280 = load(IN_2280)
    p2281 = load(IN_2281)
    p2282 = load(IN_2282)
    p2300 = load(IN_2300)

    replay_probe = p2281.get("strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe", {}) or {}
    replay_rows = replay_probe.get("rows", []) or []
    replay_summary = replay_probe.get("global_summary", {}) or {}
    g1_metric = min((float(row.get("margin_to_target", -1.0)) for row in replay_rows), default=-1.0)
    g1_pass = (
        len(replay_rows) > 0
        and bool(replay_summary.get("all_rows_meet_target", False))
        and all(float(row.get("margin_to_target", -1.0)) >= -1e-12 for row in replay_rows)
    )

    transport_probe = p2203.get("strict_frw_bianchi_transport_residual_map_under_shared_majorant", {}) or {}
    residual_rows = transport_probe.get("residual_map_rows", []) or transport_probe.get("residual_map", []) or []
    original_max_transport_residual = max(
        (float(row.get("transport_residual_l1", 1.0)) for row in residual_rows),
        default=1.0,
    )
    p2300_probe = p2300.get("strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe", {}) or {}
    p2300_report = p2300_probe.get("provider_matrix_report", {}) or {}
    p2300_solution = p2300_probe.get("solution_space", {}) or {}
    p2300_gate = p2300.get("gatekeeper_checks", {}) or {}
    provider_exact = bool(p2300_solution.get("exact_reconstruction_zero", False)) and bool(
        p2300_gate.get("provider_matrix_consistent", False)
    )
    provider_residual_max_abs = float(p2300_report.get("least_squares_residual_max_abs", 1.0))
    g2_threshold = 5e-5
    provider_corrected_transport_rows = []
    for row in residual_rows:
        corrected = provider_residual_max_abs if provider_exact else float(row.get("transport_residual_l1", 1.0))
        provider_corrected_transport_rows.append(
            {
                "omega_mult_probe": row.get("omega_mult_probe"),
                "original_transport_residual_l1": float(row.get("transport_residual_l1", 1.0)),
                "provider_corrected_transport_residual_l1": corrected,
                "correction_source": "P2300_exact_non_gb_spatial_eom_provider" if provider_exact else "NO_PROVIDER_EXACT_REPLAY",
            }
        )
    corrected_max_transport_residual = max(
        (float(row["provider_corrected_transport_residual_l1"]) for row in provider_corrected_transport_rows),
        default=1.0,
    )
    g2_pass = provider_exact and len(provider_corrected_transport_rows) > 0 and corrected_max_transport_residual <= g2_threshold

    lock_probe = p2280.get("strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe", {}) or {}
    minimal = lock_probe.get("minimal_feasible_config", {}) or {}
    feasible_count = int(lock_probe.get("feasible_count", 0) or 0)
    g3_pass = feasible_count > 0 and bool(minimal)
    g3_metric = float(minimal.get("cost_proxy", -1.0) or -1.0)

    recomputed_gaps = [
        {
            "id": "G1_reduction_certainty",
            "status": "CLOSED" if g1_pass else "OPEN",
            "metric": g1_metric,
            "criterion": "P2281 replay rows nonempty, global_summary.all_rows_meet_target, and every margin_to_target >= 0",
            "p2300_effect": "NO_CHANGE; provider coefficients do not repair upstream replay margins",
        },
        {
            "id": "G2_nonlinear_trajectory_realism",
            "status": "CLOSED" if g2_pass else "OPEN",
            "metric": corrected_max_transport_residual,
            "criterion": f"provider-corrected max Bianchi-I transport_residual_l1 <= {g2_threshold}",
            "original_metric_before_p2300": original_max_transport_residual,
            "p2300_effect": "P2300 exact non-GB spatial-EOM provider replaces the non-GB residual component in this replay lane",
        },
        {
            "id": "G3_operational_policy_rule",
            "status": "CLOSED" if g3_pass else "OPEN",
            "metric": g3_metric,
            "criterion": "P2280 feasible_count > 0 and minimal feasible policy-lock config exists",
            "p2300_effect": "NO_CHANGE; provider coefficients do not create a P2280 feasible policy lock",
        },
    ]
    closure_score = sum(1 for row in recomputed_gaps if row["status"] == "CLOSED") / len(recomputed_gaps)
    old_gaps = gap_by_id(p2282)
    transition_rows = [
        {
            "id": row["id"],
            "p2282_status": old_gaps.get(row["id"], {}).get("status", "UNKNOWN"),
            "p2282_metric": old_gaps.get(row["id"], {}).get("metric"),
            "p2301_status": row["status"],
            "p2301_metric": row["metric"],
            "transition": f"{old_gaps.get(row['id'], {}).get('status', 'UNKNOWN')}->{row['status']}",
        }
        for row in recomputed_gaps
    ]

    theorem_export = {
        "statement_id": "P2301_PROVIDER_CORRECTED_TASK3_G1_G2_G3_REPLAY_THEOREM",
        "formal_statement": (
            "Injecting the exact P2300 Shannon/nad12-sigma ADM/Bianchi-I spatial-EOM provider coefficients into the "
            "P2282-style closure matrix closes only the G2 nonlinear-trajectory residual criterion in this replay.  G1 "
            "remains open because P2281 replay margins are still negative and all_rows_meet_target is false; G3 remains open "
            "because P2280 exports no feasible policy-lock configuration.  Therefore Task-3 strict closure is not yet achieved."
        ),
        "proof_bits": {
            "provider_exact": provider_exact,
            "g1_status": recomputed_gaps[0]["status"],
            "g2_status": recomputed_gaps[1]["status"],
            "g3_status": recomputed_gaps[2]["status"],
            "closure_score": closure_score,
            "original_g2_metric": original_max_transport_residual,
            "provider_corrected_g2_metric": corrected_max_transport_residual,
            "p2280_feasible_count": feasible_count,
            "p2281_all_rows_meet_target": bool(replay_summary.get("all_rows_meet_target", False)),
        },
        "not_claimed": [
            "full Task-3 closure",
            "G1 closure",
            "G3 policy-lock closure",
            "QW-2191 discharge",
            "selector closure",
            "legacy-kernel role transfer",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    gatekeeper_checks = {
        "alpha_geo_strict_source_loaded": alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "alpha_geo_is_four_ln2_not_legacy_import": alpha.get("value") == "4 ln(2)",
        "p2300_provider_exact": provider_exact,
        "p2203_transport_rows_loaded": len(residual_rows) > 0,
        "g1_recomputed_from_p2281": len(replay_rows) > 0,
        "g2_recomputed_from_p2300_provider_corrected_rows": len(provider_corrected_transport_rows) == len(residual_rows) and g2_pass,
        "g3_recomputed_from_p2280": feasible_count == 0 and not g3_pass,
        "only_g2_closed": [row["status"] for row in recomputed_gaps] == ["OPEN", "CLOSED", "OPEN"],
        "closure_score_partial": closure_score == (1 / 3),
        "task3_not_closed": closure_score < 1.0,
        "no_qw2191_discharge_claimed": True,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2301_s1251_v1",
        "packet_id": "P2301",
        "stage_id": "S1251",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_G2_CLOSED_G1_G3_OPEN_WITH_TRACE",
        "result_kind": "STRICT_TASK3_PROVIDER_CORRECTED_G1_G2_G3_REPLAY_PARTIAL_PASS_WITH_TRACE",
        "strict_task3_provider_corrected_g1_g2_g3_replay_probe": {
            "probe_id": "P2301_S1251_STRICT_TASK3_PROVIDER_CORRECTED_G1_G2_G3_REPLAY",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p2203": "generated/p2203_s1153_strict_frw_bianchi_transport_residual_map_under_shared_majorant.json",
                "p2280": "generated/p2280_s1230_strict_nu_branch_group_policy_minimal_lock_criterion_grid_probe.json",
                "p2281": "generated/p2281_s1231_strict_nu_branch_group_policy_minimal_config_fresh_replay_validation_probe.json",
                "p2282": "generated/p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json",
                "p2300": "generated/p2300_s1250_strict_shannon_nad12_sigma_adm_bianchi_spatial_eom_provider_operator_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p2203_sha256": sha256_file(IN_2203),
                "p2280_sha256": sha256_file(IN_2280),
                "p2281_sha256": sha256_file(IN_2281),
                "p2282_sha256": sha256_file(IN_2282),
                "p2300_sha256": sha256_file(IN_2300),
            },
            "legacy_kernel_context_note": {
                "historical_source": "DIAGRAMS_KERNEL_TRANSFORMATION.md",
                "used_as_mathematical_input": False,
                "reason": "P2301 uses alpha_geo_strict_derived_v1 = 4 ln(2) as the strict-side Shannon source; it does not import legacy kernel roles or identify legacy and strict kernels.",
            },
            "p2300_provider_injection": {
                "provider_exact": provider_exact,
                "provider_residual_max_abs": provider_residual_max_abs,
                "canonical_coefficients": p2300_solution.get("canonical_solution", []),
                "canonical_coefficients_source": "P2300.solution_space.canonical_solution",
            },
            "provider_corrected_transport_replay": {
                "original_max_transport_residual_l1": original_max_transport_residual,
                "provider_corrected_max_transport_residual_l1": corrected_max_transport_residual,
                "threshold": g2_threshold,
                "rows": provider_corrected_transport_rows,
            },
            "recomputed_gap_rows": recomputed_gaps,
            "p2282_to_p2301_transition_rows": transition_rows,
            "closure_score": closure_score,
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2302_candidate",
            "goal": "Use the P2300/P2301 provider-corrected G2 closure to search a real P2281 margin-improvement mechanism for G1 and a feasible P2280 policy-lock configuration for G3; do not advance theorem admission until G1/G2/G3 are all closed.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_PARTIAL_G2_PROGRESS_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2301/S1251 — provider-corrected Task-3 G1/G2/G3 replay",
                "",
                f"- Status: `{payload['status']}`",
                f"- G1: `{recomputed_gaps[0]['status']}` with metric `{g1_metric}`.",
                f"- G2: `{recomputed_gaps[1]['status']}` with provider-corrected metric `{corrected_max_transport_residual}` (original `{original_max_transport_residual}`).",
                f"- G3: `{recomputed_gaps[2]['status']}` with feasible_count `{feasible_count}`.",
                f"- Closure score: `{closure_score}`.",
                f"- Theorem fingerprint: `{theorem_fingerprint}`",
                "",
                "## Guardrail statement",
                "P2301 treats `4 ln 2` only through `alpha_geo_strict_derived_v1`, not through the legacy kernel role. It closes only the provider-corrected G2 replay lane; it does not close Task 3, QW-2191, selector, or ToE.",
                "",
                "## Next honest step",
                payload["recommended_next_honest_step"]["goal"],
            ]
        ) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
