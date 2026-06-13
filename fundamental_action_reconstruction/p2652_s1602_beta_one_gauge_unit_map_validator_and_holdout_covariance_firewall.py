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
OUT = GEN / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.json"
MD = GEN / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.md"

P2646 = GEN / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.json"
P2647 = GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json"
P2648 = GEN / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.json"
P2651 = GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

ETA = 9.0 / 5.0
RAW_FIXTURE_BETA = 0.09270338861541028
TOL = 1e-10

REQUIRED_GAUGE_CONTRACT_KEYS = {
    "beta_gauge_status",
    "eta",
    "beta_before_gauge",
    "beta_after_gauge",
    "gauge_scale_a",
    "distance_map_direction",
    "unit_map_source",
    "unit_map_source_precommitted_before_holdout",
    "raw_distances_are_not_interpreted_as_beta_one_distances",
}

REQUIRED_UNIT_MEASUREMENT_KEYS = {
    "near",
    "far",
    "raw_near",
    "raw_far",
    "near_in_beta_one_gauge",
    "far_in_beta_one_gauge",
    "measured_tail_ratio",
    "measured_log_tail_slope",
    "ratio_standard_error",
    "slope_standard_error",
}

FORBIDDEN_UNIT_MAP_SOURCES = {
    "tail_ratio_fit",
    "holdout_fit",
    "strict_kernel_target_fit",
    "post_unblind_refit",
    "per_pair_fit",
}

NEGATIVE_EXPORT_FLAGS = [
    "real_blind_holdout_payload_loaded",
    "empirical_confirmation_claimed",
    "unit_map_source_theorem_exported",
    "target_fit_unit_map_allowed",
    "raw_distance_beta_one_substitution_allowed",
    "beta_source_exported",
    "bridge_completion_exported",
    "role_transfer_revalidated",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


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
    return {"count": len(lines), "samples": lines[:70]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "gauge_unit_map_content": (
            "gauge-fixed|gauge fixed|unit map|distance map|distance coordinates|beta=1 gauge|"
            "raw distances|gauge_scale|normalization representative"
        ),
        "holdout_covariance_content": (
            "blind holdout|holdout payload|measured_tail_ratio|measured_log_tail_slope|"
            "tail ratios.*gauge|distance coordinates are scaled|no retuning"
        ),
        "fake_unit_pass_firewall_content": (
            "fake pass|holdout fit|per-pair fit|post-unblind|raw distance|target fit|"
            "strict_kernel_target|control baseline"
        ),
        "source_nonclosure_content": (
            "target-independent beta|beta source|source theorem|role-transfer|role-bearing L_total|"
            "QW-2191|ToE closure|bridge completion"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for beta=1 gauge/unit-map holdout covariance, not packet-name search",
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


def envelope(beta: float, d: float) -> float:
    return 1.0 / (1.0 + beta * d**ETA)


def gauge_scale_to_beta_one(beta: float) -> float:
    return beta ** (1.0 / ETA)


def ratio_and_slope(beta: float, near: float, far: float) -> tuple[float, float]:
    ratio = envelope(beta, far) / envelope(beta, near)
    slope = -math.log(ratio) / math.log(far / near)
    return ratio, slope


def base_payload(contract: dict[str, Any], measurements: list[dict[str, Any]], *, unit_source: str = "external_precommitted_unit_map") -> dict[str, Any]:
    return {
        "p2646_payload_sha256": load_json(P2646).get("payload_sha256"),
        "dataset_blind_id_hash": "synthetic-fixture-not-data",
        "preprocessing_commit_hash": "synthetic-fixture-not-data",
        "phase_mask_policy_hash": "synthetic-fixture-not-data",
        "beta_gauge_contract": {**contract, "unit_map_source": unit_source},
        "measurements": measurements,
        "control_baselines": [
            {"name": "legacy_hyperbolic", "frozen_before_holdout": True, "uses_holdout_refit": False, "passes_same_thresholds": False},
            {"name": "per_pair_target_fit", "frozen_before_holdout": False, "uses_holdout_refit": True, "passes_same_thresholds": True},
        ],
        "no_retuning_attestation": True,
    }


def covariant_synthetic_payload(rows: list[dict[str, Any]], beta: float = RAW_FIXTURE_BETA) -> dict[str, Any]:
    scale = gauge_scale_to_beta_one(beta)
    measurements = []
    for row in rows:
        raw_near = row["near"] / scale
        raw_far = row["far"] / scale
        ratio, slope = ratio_and_slope(beta, raw_near, raw_far)
        measurements.append({
            "near": row["near"],
            "far": row["far"],
            "raw_near": raw_near,
            "raw_far": raw_far,
            "near_in_beta_one_gauge": scale * raw_near,
            "far_in_beta_one_gauge": scale * raw_far,
            "measured_tail_ratio": ratio,
            "measured_log_tail_slope": slope,
            "ratio_standard_error": 0.0,
            "slope_standard_error": 0.0,
        })
    contract = {
        "beta_gauge_status": "beta=1 gauge-fixed working normalization",
        "eta": ETA,
        "beta_before_gauge": beta,
        "beta_after_gauge": 1.0,
        "gauge_scale_a": scale,
        "distance_map_direction": "d_beta_one = gauge_scale_a * d_raw",
        "unit_map_source_precommitted_before_holdout": True,
        "raw_distances_are_not_interpreted_as_beta_one_distances": True,
    }
    return base_payload(contract, measurements)


def raw_substitution_fake_payload(rows: list[dict[str, Any]], beta: float = RAW_FIXTURE_BETA) -> dict[str, Any]:
    scale = gauge_scale_to_beta_one(beta)
    measurements = []
    for row in rows:
        ratio, slope = ratio_and_slope(1.0, row["near"], row["far"])
        measurements.append({
            "near": row["near"],
            "far": row["far"],
            "raw_near": row["near"],
            "raw_far": row["far"],
            "near_in_beta_one_gauge": row["near"],
            "far_in_beta_one_gauge": row["far"],
            "measured_tail_ratio": ratio,
            "measured_log_tail_slope": slope,
            "ratio_standard_error": 0.0,
            "slope_standard_error": 0.0,
        })
    contract = {
        "beta_gauge_status": "beta=1 asserted on raw distances",
        "eta": ETA,
        "beta_before_gauge": beta,
        "beta_after_gauge": 1.0,
        "gauge_scale_a": scale,
        "distance_map_direction": "d_beta_one = gauge_scale_a * d_raw",
        "unit_map_source_precommitted_before_holdout": True,
        "raw_distances_are_not_interpreted_as_beta_one_distances": False,
    }
    return base_payload(contract, measurements)


def target_fit_fake_payload(rows: list[dict[str, Any]]) -> dict[str, Any]:
    payload = covariant_synthetic_payload(rows)
    payload["beta_gauge_contract"]["unit_map_source_precommitted_before_holdout"] = False
    payload["beta_gauge_contract"]["raw_distances_are_not_interpreted_as_beta_one_distances"] = True
    payload["beta_gauge_contract"]["unit_map_source"] = "tail_ratio_fit"
    return payload


def validate_unit_map_payload(candidate: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    contract = candidate.get("beta_gauge_contract", {}) if isinstance(candidate.get("beta_gauge_contract", {}), dict) else {}
    missing_contract = sorted(REQUIRED_GAUGE_CONTRACT_KEYS - set(contract))
    measurements = candidate.get("measurements", []) if isinstance(candidate.get("measurements", []), list) else []
    expected_pairs = {(row["near"], row["far"]) for row in rows}
    measurement_errors: list[dict[str, Any]] = []
    covariance_errors: list[dict[str, Any]] = []
    seen_pairs = set()

    beta_before = contract.get("beta_before_gauge")
    beta_after = contract.get("beta_after_gauge")
    eta = contract.get("eta")
    scale = contract.get("gauge_scale_a")
    source = contract.get("unit_map_source")
    source_forbidden = source in FORBIDDEN_UNIT_MAP_SOURCES
    source_precommitted = contract.get("unit_map_source_precommitted_before_holdout") is True
    raw_distance_guard = contract.get("raw_distances_are_not_interpreted_as_beta_one_distances") is True

    scale_equation_ok = False
    if all(isinstance(value, (int, float)) for value in [beta_before, beta_after, eta, scale]) and beta_after:
        scale_equation_ok = abs(float(scale) ** float(eta) - float(beta_before) / float(beta_after)) < TOL

    for index, measurement in enumerate(measurements):
        missing = sorted(REQUIRED_UNIT_MEASUREMENT_KEYS - set(measurement))
        pair = (measurement.get("near"), measurement.get("far"))
        seen_pairs.add(pair)
        pair_errors = []
        if pair not in expected_pairs:
            pair_errors.append("pair_not_preregistered")
        if missing:
            pair_errors.append("missing_required_measurement_keys")
        if isinstance(scale, (int, float)) and not missing:
            near_cov = abs(float(scale) * float(measurement["raw_near"]) - float(measurement["near_in_beta_one_gauge"])) < TOL
            far_cov = abs(float(scale) * float(measurement["raw_far"]) - float(measurement["far_in_beta_one_gauge"])) < TOL
            near_locked = abs(float(measurement["near_in_beta_one_gauge"]) - float(measurement["near"])) < TOL
            far_locked = abs(float(measurement["far_in_beta_one_gauge"]) - float(measurement["far"])) < TOL
            if not (near_cov and far_cov and near_locked and far_locked):
                covariance_errors.append({
                    "index": index,
                    "pair": pair,
                    "near_covariant": near_cov,
                    "far_covariant": far_cov,
                    "near_matches_locked_p2646_pair": near_locked,
                    "far_matches_locked_p2646_pair": far_locked,
                })
        if pair_errors:
            measurement_errors.append({"index": index, "pair": pair, "missing_keys": missing, "errors": pair_errors})

    all_pairs_present = expected_pairs == seen_pairs
    schema_ready = not missing_contract and not measurement_errors and all_pairs_present
    unit_map_admissible = schema_ready and scale_equation_ok and not covariance_errors and source_precommitted and not source_forbidden and raw_distance_guard
    return {
        "schema_ready": schema_ready,
        "missing_gauge_contract_keys": missing_contract,
        "measurement_errors": measurement_errors,
        "all_preregistered_pairs_present": all_pairs_present,
        "scale_equation_ok": scale_equation_ok,
        "covariance_errors": covariance_errors,
        "unit_map_source": source,
        "unit_map_source_precommitted_before_holdout": source_precommitted,
        "unit_map_source_forbidden_as_target_fit": source_forbidden,
        "raw_distance_guard_declared": raw_distance_guard,
        "unit_map_admissible_for_p2647_p2648_execution": unit_map_admissible,
        "real_payload_loaded": False,
    }


def firewall(rows: list[dict[str, Any]]) -> dict[str, Any]:
    covariant = validate_unit_map_payload(covariant_synthetic_payload(rows), rows)
    raw_fake = validate_unit_map_payload(raw_substitution_fake_payload(rows), rows)
    target_fake = validate_unit_map_payload(target_fit_fake_payload(rows), rows)
    return {
        "covariant_synthetic_fixture": covariant,
        "raw_distance_beta_one_substitution_fixture": raw_fake,
        "target_fit_unit_map_fixture": target_fake,
        "firewall_passes": (
            covariant["unit_map_admissible_for_p2647_p2648_execution"]
            and not raw_fake["unit_map_admissible_for_p2647_p2648_execution"]
            and not target_fake["unit_map_admissible_for_p2647_p2648_execution"]
        ),
        "blocked_fake_passes": [
            "setting beta=1 while leaving raw distances unscaled",
            "deriving the unit map from holdout tail-ratio/strict-target fit",
            "per-pair post-unblind unit-map retuning",
        ],
    }


def upstream_consistency() -> dict[str, Any]:
    p2647 = load_json(P2647)
    p2648 = load_json(P2648)
    p2651 = load_json(P2651)
    return {
        "p2647_schema_ready_but_no_real_holdout": p2647.get("closure_decision", {}).get("real_blind_holdout_executed") is False,
        "p2648_margin_ready_but_no_real_holdout": p2648.get("closure_decision", {}).get("real_blind_holdout_executed") is False,
        "p2651_allows_beta_one_only_as_working_gauge": p2651.get("closure_decision", {}).get("beta_one_working_gauge_allowed") is True,
        "p2651_blocks_beta_source": p2651.get("closure_decision", {}).get("beta_source_exported_now") is False,
        "p2651_requires_unit_map_bookkeeping": "unit-map bookkeeping" in json.dumps(p2651, sort_keys=True),
    }


def closure_decision(fw: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "GAUGE_UNIT_MAP_VALIDATOR_READY_BUT_NO_REAL_PAYLOAD_NO_UNIT_SOURCE_THEOREM_NO_CLOSURE",
        "professorial_verdict": (
            "P2652 turns the P2651 beta=1 gauge contract into an executable covariance firewall for P2647/P2648.  "
            "A holdout payload may be evaluated in the beta=1 working gauge only if it carries a precommitted distance/unit map whose scaled distances match the locked P2646 pairs.  "
            "The covariant synthetic fixture passes that bookkeeping check, while raw-distance beta=1 substitution and target/holdout-fitted unit maps are rejected.  "
            "This is still validator readiness only: it supplies no real data and no typed nadsoliton unit-source theorem."
        ),
        "professorial_closure_path": [
            "Require every real P2647/P2648 payload to include beta_gauge_contract and per-measurement raw-to-beta-one distance coordinates.",
            "Reject any payload whose unit map is fitted from the holdout tail ratios, from the strict target, or per pair after unblinding.",
            "Run the statistical margin rule only after the gauge/unit covariance validator passes.",
            "Keep beta=1 labelled as gauge-fixed working normalization unless a typed nadsoliton metric/UV source theorem is proved independently.",
            "If real data pass P2647/P2648/P2652, rerun only the modified/compressed successor role row; do not revive unchanged legacy inverse-hierarchy or L_total closure.",
        ],
        "next_honest_step": (
            "Prepare a real blinded payload with a precommitted raw-to-beta=1 unit map, run P2652 first, then P2647/P2648 exactly once; "
            "in parallel, attempt the independent typed nadsoliton metric/UV source theorem rather than fitting the unit map from the compression target."
        ),
        "unit_map_firewall_passes": fw["firewall_passes"],
        "real_blind_holdout_payload_loaded": False,
        "unit_map_source_theorem_exported_now": False,
        "empirical_confirmation_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    fw = payload["gauge_unit_map_firewall"]
    lines = [
        "# P2652/S1602 beta=1 gauge unit-map validator and holdout covariance firewall",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps gauge/unit-map, holdout covariance, fake unit-pass firewall, and source nonclosure content before adding the validator.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Gauge/unit-map validator schema",
        "",
        f"Required gauge contract keys: `{sorted(REQUIRED_GAUGE_CONTRACT_KEYS)}`.",
        f"Required unit-aware measurement keys: `{sorted(REQUIRED_UNIT_MEASUREMENT_KEYS)}`.",
        f"Forbidden unit-map sources: `{sorted(FORBIDDEN_UNIT_MAP_SOURCES)}`.",
        "",
        "## Firewall result",
        "",
        f"Covariant synthetic fixture admissible? `{fw['covariant_synthetic_fixture']['unit_map_admissible_for_p2647_p2648_execution']}`.",
        f"Raw-distance beta=1 substitution admissible? `{fw['raw_distance_beta_one_substitution_fixture']['unit_map_admissible_for_p2647_p2648_execution']}`.",
        f"Target-fit unit map admissible? `{fw['target_fit_unit_map_fixture']['unit_map_admissible_for_p2647_p2648_execution']}`.",
        f"Firewall passes? `{fw['firewall_passes']}`.",
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Real blind holdout payload loaded? `{decision['real_blind_holdout_payload_loaded']}`.",
        f"Unit-map source theorem exported now? `{decision['unit_map_source_theorem_exported_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
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
    rows = p2646_rows()
    fw = firewall(rows)
    decision = closure_decision(fw)
    payload: dict[str, Any] = {
        "status": "P2652_GAUGE_UNIT_MAP_VALIDATOR_READY_NO_REAL_PAYLOAD_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2646_PREREGISTRATION": sha256_file(P2646),
            "P2647_HOLDOUT_HARNESS": sha256_file(P2647),
            "P2648_MARGIN_RULE": sha256_file(P2648),
            "P2651_GAUGE_CONTRACT": sha256_file(P2651),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "required_gauge_unit_payload_schema": {
            "gauge_contract_keys": sorted(REQUIRED_GAUGE_CONTRACT_KEYS),
            "unit_aware_measurement_keys": sorted(REQUIRED_UNIT_MEASUREMENT_KEYS),
            "forbidden_unit_map_sources": sorted(FORBIDDEN_UNIT_MAP_SOURCES),
        },
        "gauge_unit_map_firewall": fw,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2652/S1602 beta=1 gauge unit-map validator guard",
        "## P2652/S1602 beta=1 gauge unit-map validator guard\n\n"
        "`P2652/S1602` turns the `P2651/S1601` beta=1 working-gauge contract into a unit-map/covariance validator for `P2647/P2648` payloads.  "
        "A holdout payload must provide a precommitted raw-to-beta=1 distance map and cannot obtain that map from holdout tail-ratio fitting, strict-target fitting, or post-unblind per-pair retuning.  "
        "This exports validator readiness only: no real empirical confirmation, no typed unit-source theorem, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2652/S1602 beta=1 gauge unit-map Ltotal guard",
        "## P2652/S1602 beta=1 gauge unit-map Ltotal guard\n\n"
        "`P2652/S1602` does not re-enable `L_total`: it blocks raw-distance beta=1 substitution and target-fitted unit maps before any frozen holdout is interpreted, but the variational damping term still needs a typed nadsoliton metric/UV source theorem, real blind data, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
