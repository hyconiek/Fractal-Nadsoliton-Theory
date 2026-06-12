#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json"
MD = GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.md"

P2646 = GEN / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "blind_holdout_data_loaded",
    "empirical_confirmation_claimed",
    "synthetic_fixture_promoted_to_evidence",
    "kernel_retuning_allowed",
    "control_baseline_defeated_by_data",
    "beta_source_exported",
    "role_transfer_revalidated",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]

REQUIRED_BLIND_SCHEMA_KEYS = {
    "p2646_payload_sha256",
    "dataset_blind_id_hash",
    "preprocessing_commit_hash",
    "phase_mask_policy_hash",
    "measurements",
    "control_baselines",
    "no_retuning_attestation",
}

REQUIRED_MEASUREMENT_KEYS = {
    "near",
    "far",
    "measured_tail_ratio",
    "measured_log_tail_slope",
    "ratio_standard_error",
    "slope_standard_error",
    "n_effective_bins",
}

REQUIRED_CONTROL_KEYS = {
    "name",
    "frozen_before_holdout",
    "uses_holdout_refit",
    "passes_same_thresholds",
}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean", "-g", "*.csv", "-g", "*.tsv",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "blind_holdout_execution_content": (
            "blind holdout execution|holdout execution|dataset_blind_id|blind_id|no-retuning|no retuning|"
            "preprocessing_commit|phase mask policy|phase_mask"
        ),
        "tail_ratio_measurement_schema_content": (
            "measured_tail_ratio|denominator/envelope tail ratio|tail-ratio/log-slope|log-tail slope|"
            "n_effective_bins|ratio_standard_error|slope_standard_error"
        ),
        "fake_pass_or_control_firewall_content": (
            "fake pass|synthetic fixture|control baseline|matched exponential|matched spline|holdout refit|"
            "adversarial controls|kernel retuning"
        ),
        "source_nonclosure_content": (
            "beta source|target-independent|role-transfer|role-bearing L_total|QW-2191|ToE closure|"
            "empirical confirmation"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for blind-holdout execution harness and fake-pass firewall, not packet-name search",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(["git", "log", "-n", "12", "--oneline", "--name-only"], cwd=REPO, text=True, stdout=subprocess.PIPE, check=True)
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def p2646_rows() -> list[dict[str, Any]]:
    payload = load_json(P2646)
    return payload.get("preregistered_frozen_kernel_tests", {}).get("rows", [])


def evaluate_measurement(row: dict[str, Any], measurement: dict[str, Any]) -> dict[str, Any]:
    ratio_margin = row["midpoint_ratio_threshold"] - measurement["measured_tail_ratio"]
    slope_margin = measurement["measured_log_tail_slope"] - row["midpoint_slope_threshold"]
    return {
        "near": row["near"],
        "far": row["far"],
        "ratio_pass": ratio_margin > 0,
        "slope_pass": slope_margin > 0,
        "ratio_margin_to_locked_threshold": ratio_margin,
        "slope_margin_to_locked_threshold": slope_margin,
        "combined_pass": ratio_margin > 0 and slope_margin > 0,
    }


def measurement_from_row(row: dict[str, Any], model: str) -> dict[str, Any]:
    if model == "strict_fixture":
        ratio = row["strict_denominator_tail_ratio"]
        slope = row["strict_log_tail_slope"]
    elif model == "legacy_fixture":
        ratio = row["legacy_hyperbolic_tail_ratio"]
        slope = row["legacy_log_tail_slope"]
    elif model == "midpoint_fixture":
        ratio = row["midpoint_ratio_threshold"]
        slope = row["midpoint_slope_threshold"]
    else:
        raise ValueError(model)
    return {
        "near": row["near"],
        "far": row["far"],
        "measured_tail_ratio": ratio,
        "measured_log_tail_slope": slope,
        "ratio_standard_error": 0.0,
        "slope_standard_error": 0.0,
        "n_effective_bins": math.inf,
    }


def evaluate_fixture(model: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    pair_results = [evaluate_measurement(row, measurement_from_row(row, model)) for row in rows]
    return {
        "model": model,
        "pair_results": pair_results,
        "all_pairs_pass": all(result["combined_pass"] for result in pair_results),
        "minimum_ratio_margin": min(result["ratio_margin_to_locked_threshold"] for result in pair_results),
        "minimum_slope_margin": min(result["slope_margin_to_locked_threshold"] for result in pair_results),
    }


def validate_blind_payload(candidate: dict[str, Any], p2646_hash: str | None, rows: list[dict[str, Any]]) -> dict[str, Any]:
    missing_top = sorted(REQUIRED_BLIND_SCHEMA_KEYS - set(candidate))
    measurement_errors: list[dict[str, Any]] = []
    measurements = candidate.get("measurements", []) if isinstance(candidate.get("measurements", []), list) else []
    expected_pairs = {(row["near"], row["far"]) for row in rows}
    seen_pairs = set()
    for index, measurement in enumerate(measurements):
        missing = sorted(REQUIRED_MEASUREMENT_KEYS - set(measurement))
        pair = (measurement.get("near"), measurement.get("far"))
        seen_pairs.add(pair)
        if missing or pair not in expected_pairs:
            measurement_errors.append({"index": index, "pair": pair, "missing_keys": missing, "pair_is_preregistered": pair in expected_pairs})
    control_errors: list[dict[str, Any]] = []
    controls = candidate.get("control_baselines", []) if isinstance(candidate.get("control_baselines", []), list) else []
    for index, control in enumerate(controls):
        missing = sorted(REQUIRED_CONTROL_KEYS - set(control))
        if missing:
            control_errors.append({"index": index, "missing_keys": missing})
    hash_matches = candidate.get("p2646_payload_sha256") == p2646_hash
    no_retuning = candidate.get("no_retuning_attestation") is True
    all_pairs_present = expected_pairs == seen_pairs
    controls_frozen = bool(controls) and all(control.get("frozen_before_holdout") is True and control.get("uses_holdout_refit") is False for control in controls)
    controls_do_not_remove_advantage = bool(controls) and all(control.get("passes_same_thresholds") is False for control in controls)
    schema_ready = not missing_top and not measurement_errors and not control_errors
    admissible_for_execution = schema_ready and hash_matches and no_retuning and all_pairs_present and controls_frozen and controls_do_not_remove_advantage
    return {
        "schema_ready": schema_ready,
        "missing_top_level_keys": missing_top,
        "measurement_errors": measurement_errors,
        "control_errors": control_errors,
        "p2646_hash_matches": hash_matches,
        "no_retuning_attested": no_retuning,
        "all_preregistered_pairs_present": all_pairs_present,
        "controls_frozen_without_holdout_refit": controls_frozen,
        "controls_do_not_remove_frozen_kernel_advantage": controls_do_not_remove_advantage,
        "admissible_for_execution": admissible_for_execution,
    }


def fake_pass_firewall(rows: list[dict[str, Any]]) -> dict[str, Any]:
    strict = evaluate_fixture("strict_fixture", rows)
    legacy = evaluate_fixture("legacy_fixture", rows)
    midpoint = evaluate_fixture("midpoint_fixture", rows)
    return {
        "purpose": "unit-test the locked discriminator without promoting synthetic fixtures to empirical evidence",
        "strict_fixture_expected_pass": strict,
        "legacy_fixture_expected_fail": legacy,
        "midpoint_fixture_expected_fail_by_strict_inequality": midpoint,
        "free_per_pair_fit_warning": (
            "A per-pair spline/power-law/exponential fitted on the holdout has at least one effective degree of freedom per audited pair and can mimic thresholds; "
            "therefore such a control is not admissible evidence and is explicitly blocked by the schema."
        ),
        "firewall_passes": strict["all_pairs_pass"] and not legacy["all_pairs_pass"] and not midpoint["all_pairs_pass"],
    }


def blind_data_presence_audit() -> dict[str, Any]:
    candidate = rg_count(
        "dataset_blind_id_hash|measured_tail_ratio|phase_mask_policy_hash|no_retuning_attestation|"
        "ratio_standard_error|slope_standard_error"
    )
    generated_hits = subprocess.run(
        ["rg", "-n", "dataset_blind_id_hash|measured_tail_ratio|phase_mask_policy_hash|no_retuning_attestation", "fundamental_action_reconstruction/generated", "-g", "*.json", "-g", "*.md"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return {
        "non_generated_candidate_hits": candidate,
        "generated_protocol_hits_count": len([line for line in generated_hits.stdout.splitlines() if line]),
        "real_blind_holdout_payload_found": False,
        "reason": "No external blinded measurement payload is supplied to this script; generated protocol text and synthetic fixtures are not data.",
    }


def closure_decision(firewall: dict[str, Any], data_audit: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "EXECUTION_HARNESS_READY_BUT_BLIND_HOLDOUT_NOT_EXECUTED__NO_EMPIRICAL_CONFIRMATION_NO_SOURCE_EXPORT",
        "professorial_verdict": (
            "P2647 is the honest post-preregistration move: it turns P2646 from a static threshold table into an executable blind-holdout harness, "
            "validates the required schema, and proves by synthetic controls that strict frozen predictions would pass while legacy and exact-midpoint fixtures fail. "
            "Because no real blinded measurement payload is loaded, this remains readiness/fake-pass-firewall evidence only, not empirical confirmation."
        ),
        "professorial_closure_path": [
            "Obtain or export a real blinded denominator/envelope measurement payload with all P2647 schema keys and the exact P2646 payload hash.",
            "Run the harness once on that payload with no beta/eta/phase retuning and with frozen controls recorded before unblinding.",
            "If any control baseline fitted without holdout refit removes the frozen-kernel advantage, mark the empirical compression signature failed, not merely inconclusive.",
            "If the holdout passes, rerun only the modified/compressed successor row of the role-transfer matrix; do not resurrect unchanged inverse-hierarchy.",
            "Keep beta-source and QW-2191/source work independent: empirical success can prioritize those proofs but cannot replace them.",
        ],
        "next_honest_step": (
            "Run P2647 on a real blinded payload matching the locked schema, or build the target-independent beta-source theorem. "
            "Do not treat the synthetic strict fixture as data and do not reopen L_total until blind data, controls, beta source, and selector/source blockers are all discharged."
        ),
        "harness_firewall_passes": firewall["firewall_passes"],
        "real_blind_holdout_executed": data_audit["real_blind_holdout_payload_found"],
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    firewall = payload["fake_pass_firewall"]
    decision = payload["closure_decision"]
    lines = [
        "# P2647/S1597 frozen-kernel blind-holdout execution harness and fake-pass firewall",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps for blind holdout execution, measured tail-ratio/log-slope schema, fake-pass/control-baseline language, and source nonclosure guardrails before adding the harness.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Harness result",
        "",
        f"P2646 payload hash: `{payload['p2646_payload_sha256']}`.",
        f"Real blind holdout executed? `{decision['real_blind_holdout_executed']}`.",
        f"Fake-pass firewall passes? `{decision['harness_firewall_passes']}`.",
        "",
        "| fixture | all pairs pass | min ratio margin | min slope margin |",
        "| --- | ---: | ---: | ---: |",
    ])
    for key in ["strict_fixture_expected_pass", "legacy_fixture_expected_fail", "midpoint_fixture_expected_fail_by_strict_inequality"]:
        result = firewall[key]
        lines.append(
            f"| `{result['model']}` | `{result['all_pairs_pass']}` | "
            f"`{result['minimum_ratio_margin']:.12f}` | `{result['minimum_slope_margin']:.12f}` |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Full kernel now? `{decision['full_kernel_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Professorial closure path",
        "",
    ])
    for item in decision["professorial_closure_path"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Next honest step",
        "",
        decision["next_honest_step"],
        "",
        "## Negative exports",
        "",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2646_payload = load_json(P2646)
    rows = p2646_rows()
    p2646_hash = p2646_payload.get("payload_sha256")
    firewall = fake_pass_firewall(rows)
    data_audit = blind_data_presence_audit()
    empty_payload_validation = validate_blind_payload({}, p2646_hash, rows)
    decision = closure_decision(firewall, data_audit)
    payload: dict[str, Any] = {
        "status": "P2647_EXECUTION_HARNESS_READY_BUT_BLIND_HOLDOUT_NOT_EXECUTED_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2646_PREREGISTRATION": sha256_file(P2646),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "p2646_payload_sha256": p2646_hash,
        "required_blind_payload_schema": {
            "top_level_keys": sorted(REQUIRED_BLIND_SCHEMA_KEYS),
            "measurement_keys": sorted(REQUIRED_MEASUREMENT_KEYS),
            "control_baseline_keys": sorted(REQUIRED_CONTROL_KEYS),
        },
        "empty_payload_validation_witness": empty_payload_validation,
        "blind_data_presence_audit": data_audit,
        "fake_pass_firewall": firewall,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2647/S1597 frozen-kernel blind-holdout harness guard",
        "## P2647/S1597 frozen-kernel blind-holdout harness guard\n\n"
        "`P2647/S1597` converts the `P2646/S1596` frozen compression preregistration into an executable schema/harness and fake-pass firewall.  "
        "Synthetic strict fixtures pass the locked tail-ratio/log-slope inequalities while legacy and midpoint fixtures fail, but no real blinded measurement payload is loaded.  "
        "Therefore this exports harness readiness only: no empirical confirmation, no beta source, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2647/S1597 frozen-kernel blind-holdout harness Ltotal guard",
        "## P2647/S1597 frozen-kernel blind-holdout harness Ltotal guard\n\n"
        "`P2647/S1597` does not re-enable `L_total`: it supplies a blind-payload schema, discriminator evaluator, and fake-pass firewall for the modified compression successor, "
        "but a role-bearing variational term still requires actual blind data, frozen controls, a target-independent beta/source theorem, and the later role-transfer rerun.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
