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
OUT = GEN / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.json"
MD = GEN / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.md"

P2649 = GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json"
P2651 = GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json"
P2652 = GEN / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.json"
P2653 = GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

ETA = 9.0 / 5.0
REFERENCE_BETA = 1.0
AUDITED_BETAS = [0.01, 0.09270338861541028, 1.0, 1.1473957999384183, 4.0]
REFERENCE_DISTANCES = [1.0, 2.0, 7.0, 12.0]
TAIL_PAIRS = [(1.0, 7.0), (1.0, 12.0), (2.0, 8.0), (3.0, 9.0)]
TOL = 1e-12

NEGATIVE_EXPORT_FLAGS = [
    "scale_invariant_beta_selector_exported",
    "canonical_unit_exported",
    "typed_metric_uv_source_theorem_exported",
    "target_independent_beta_source_exported",
    "raw_coordinate_anchor_promoted_to_source",
    "empirical_holdout_promoted_to_source",
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
        "scale_invariant_selector_content": (
            "scale-invariant|scale invariant|normalization orbit|scale orbit|equivariant|"
            "beta selector|selector.*beta|representative"
        ),
        "typed_unit_source_content": (
            "typed nadsoliton metric|UV unit|canonical unit|distance functional|"
            "nadsoliton-selected|unit selector"
        ),
        "raw_anchor_warning_content": (
            "raw distance|raw coordinate|unit map|gauge-fixed|gauge fixed|target-fitted|"
            "holdout.*unit"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|"
            "beta source|source theorem"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for scale-invariant beta selector no-go, not packet-name search",
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


def envelope(beta: float, d: float) -> float:
    return 1.0 / (1.0 + beta * d**ETA)


def scale_to_reference(beta: float) -> float:
    return (beta / REFERENCE_BETA) ** (1.0 / ETA)


def log_slope(beta: float, near: float, far: float) -> float:
    ratio = envelope(beta, far) / envelope(beta, near)
    return -math.log(ratio) / math.log(far / near)


def invariant_feature_vector(beta: float) -> list[float]:
    scale = scale_to_reference(beta)
    distances = [d / scale for d in REFERENCE_DISTANCES]
    pairs = [(near / scale, far / scale) for near, far in TAIL_PAIRS]
    envelope_features = [envelope(beta, d) for d in distances]
    ratio_features = [envelope(beta, far) / envelope(beta, near) for near, far in pairs]
    slope_features = [log_slope(beta, near, far) for near, far in pairs]
    x_features = [beta * d**ETA for d in distances]
    return envelope_features + ratio_features + slope_features + x_features


def raw_anchor_feature_vector(beta: float) -> list[float]:
    return [envelope(beta, d) for d in REFERENCE_DISTANCES] + [log_slope(beta, near, far) for near, far in TAIL_PAIRS]


def feature_rank_no_go() -> dict[str, Any]:
    reference = invariant_feature_vector(REFERENCE_BETA)
    invariant_rows = []
    raw_rows = []
    for beta in AUDITED_BETAS:
        invariant = invariant_feature_vector(beta)
        raw = raw_anchor_feature_vector(beta)
        invariant_rows.append({
            "beta": beta,
            "scale_to_beta_one": scale_to_reference(beta),
            "max_abs_difference_from_beta_one_invariant_features": max(abs(a - b) for a, b in zip(invariant, reference)),
            "invariant_feature_vector": invariant,
        })
        raw_rows.append({
            "beta": beta,
            "max_abs_difference_from_beta_one_raw_anchor_features": max(abs(a - b) for a, b in zip(raw, raw_anchor_feature_vector(REFERENCE_BETA))),
            "raw_anchor_feature_vector": raw,
        })
    invariant_rank = 1 if all(row["max_abs_difference_from_beta_one_invariant_features"] < TOL for row in invariant_rows) else len(AUDITED_BETAS)
    raw_distinguishes = any(row["max_abs_difference_from_beta_one_raw_anchor_features"] > 1e-6 for row in raw_rows if abs(row["beta"] - 1.0) > 1e-12)
    return {
        "statement": "Any selector depending only on orbit-invariant envelope/tail/log-slope/x features is constant along the beta-distance rescaling orbit and cannot select beta=1 as a unique representative.",
        "invariant_rows": invariant_rows,
        "raw_anchor_rows": raw_rows,
        "invariant_feature_rank_on_audited_orbit": invariant_rank,
        "max_invariant_difference": max(row["max_abs_difference_from_beta_one_invariant_features"] for row in invariant_rows),
        "raw_anchor_features_distinguish_representatives": raw_distinguishes,
        "raw_anchor_warning": "Raw-coordinate features distinguish beta representatives only after an external coordinate unit has already been chosen; that is a gauge/unit premise, not an invariant beta source.",
        "scale_invariant_selector_exists_on_audited_features": False,
    }


def selector_candidate_matrix(no_go: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = [
        {
            "candidate": "orbit_invariant_envelope_grid_selector",
            "uses_only_scale_invariant_features": True,
            "uses_external_unit_anchor": False,
            "unique_beta_one_on_audited_orbit": False,
            "verdict": "blocked because the invariant feature rank along the audited orbit is one",
        },
        {
            "candidate": "orbit_invariant_tail_ratio_log_slope_selector",
            "uses_only_scale_invariant_features": True,
            "uses_external_unit_anchor": False,
            "unique_beta_one_on_audited_orbit": False,
            "verdict": "blocked because covariant tail ratios/log-slopes are identical across representatives",
        },
        {
            "candidate": "raw_coordinate_envelope_anchor_selector",
            "uses_only_scale_invariant_features": False,
            "uses_external_unit_anchor": True,
            "unique_beta_one_on_audited_orbit": no_go["raw_anchor_features_distinguish_representatives"],
            "verdict": "can distinguish representatives only by importing a raw coordinate unit, so it is a gauge premise rather than a source theorem",
        },
        {
            "candidate": "p2652_precommitted_unit_map_selector",
            "uses_only_scale_invariant_features": False,
            "uses_external_unit_anchor": True,
            "unique_beta_one_on_audited_orbit": True,
            "verdict": "validates a declared map but does not source the map from nadsoliton dynamics",
        },
        {
            "candidate": "typed_metric_uv_plus_operator_identity_selector",
            "uses_only_scale_invariant_features": False,
            "uses_external_unit_anchor": False,
            "unique_beta_one_on_audited_orbit": None,
            "verdict": "the only admissible theorem target: first source the unit from nadsoliton dynamics, then prove a dimensionless identity with unique beta=1",
        },
    ]
    for row in candidates:
        row["passes_as_target_independent_beta_source_now"] = False
    return candidates


def upstream_consistency() -> dict[str, Any]:
    p2649 = load_json(P2649)
    p2651 = load_json(P2651)
    p2652 = load_json(P2652)
    p2653 = load_json(P2653)
    return {
        "p2649_beta_routes_blocked": p2649.get("closure_decision", {}).get("passing_beta_source_routes") == [],
        "p2651_beta_one_only_working_gauge": p2651.get("closure_decision", {}).get("beta_one_working_gauge_allowed") is True,
        "p2652_unit_map_source_not_exported": p2652.get("closure_decision", {}).get("unit_map_source_theorem_exported_now") is False,
        "p2653_typed_metric_uv_source_not_exported": p2653.get("closure_decision", {}).get("typed_metric_uv_source_theorem_exported_now") is False,
        "p2653_scale_orbit_equivalence_verified": p2653.get("closure_decision", {}).get("scale_orbit_equivalence_verified") is True,
    }


def closure_decision(no_go: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_target_independent_beta_source_now"]]
    return {
        "decision": "SCALE_INVARIANT_BETA_SELECTOR_NO_GO__RAW_ANCHORS_ARE_GAUGE_PREMISES_NOT_SOURCES",
        "professorial_verdict": (
            "P2654 proves the next no-go around the P2653 obligation matrix: scale-invariant functionals of the strict denominator data are constant on the beta-distance orbit, so they cannot select beta=1.  "
            "Raw-coordinate anchors can distinguish representatives only by assuming a unit before the selector acts, which is exactly the missing typed metric/UV premise.  "
            "Therefore the proof route must source the unit from nadsoliton dynamics first; no invariant beta selector, canonical unit, or role-bearing dynamics is exported here."
        ),
        "professorial_closure_path": [
            "Stop searching for a scale-invariant beta=1 selector inside envelope/tail/log-slope data alone; the audited feature rank is one along the orbit.",
            "If a raw-coordinate anchor is used, label it as an external gauge/unit premise until a typed nadsoliton metric/UV theorem derives it.",
            "The only proof-grade route left is to construct the typed metric/UV unit first and then add a separate dimensionless operator/conservation identity with unique beta=1.",
            "Empirical P2652/P2647/P2648 remains support only and must not be reinterpreted as the missing invariant selector.",
        ],
        "next_honest_step": (
            "Attempt a constructive typed nadsoliton metric/UV theorem, not another invariant beta selector: define the nadsoliton state space and metric, derive a UV unit from its dynamics, then test a beta-selecting operator/conservation identity after the scale orbit is quotiented."
        ),
        "passing_beta_source_candidates": passing,
        "scale_invariant_selector_exists_now": no_go["scale_invariant_selector_exists_on_audited_features"],
        "raw_anchor_promoted_to_source_now": False,
        "canonical_unit_exported_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    no_go = payload["feature_rank_no_go"]
    decision = payload["closure_decision"]
    lines = [
        "# P2654/S1604 scale-invariant beta selector no-go and equivariant functional rank certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps scale-invariant selector, typed unit source, raw anchor warning, and nonclosure guard content before adding the no-go.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Feature-rank no-go",
        "",
        no_go["statement"],
        f"Invariant feature rank on audited orbit: `{no_go['invariant_feature_rank_on_audited_orbit']}`.",
        f"Max invariant feature difference from beta=1: `{no_go['max_invariant_difference']:.3e}`.",
        f"Raw-anchor features distinguish representatives? `{no_go['raw_anchor_features_distinguish_representatives']}`.",
        no_go["raw_anchor_warning"],
        "",
        "## Selector candidate matrix",
        "",
        "| candidate | invariant only? | external unit anchor? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | --- |",
    ])
    for row in payload["selector_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['uses_only_scale_invariant_features']}` | "
            f"`{row['uses_external_unit_anchor']}` | `{row['passes_as_target_independent_beta_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing beta-source candidates: `{decision['passing_beta_source_candidates']}`.",
        f"Scale-invariant selector exists now? `{decision['scale_invariant_selector_exists_now']}`.",
        f"Raw anchor promoted to source now? `{decision['raw_anchor_promoted_to_source_now']}`.",
        f"Canonical unit exported now? `{decision['canonical_unit_exported_now']}`.",
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
    no_go = feature_rank_no_go()
    matrix = selector_candidate_matrix(no_go)
    decision = closure_decision(no_go, matrix)
    payload: dict[str, Any] = {
        "status": "P2654_SCALE_INVARIANT_BETA_SELECTOR_NO_GO_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2649_BETA_ROUTE_MATRIX": sha256_file(P2649),
            "P2651_GAUGE_CONTRACT": sha256_file(P2651),
            "P2652_UNIT_MAP_VALIDATOR": sha256_file(P2652),
            "P2653_TYPED_METRIC_UV_OBLIGATION": sha256_file(P2653),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "feature_rank_no_go": no_go,
        "selector_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2654/S1604 scale-invariant beta selector no-go guard",
        "## P2654/S1604 scale-invariant beta selector no-go guard\n\n"
        "`P2654/S1604` audits whether a scale-invariant functional of strict envelope/tail/log-slope data can select `beta=1`.  The invariant feature rank on the audited beta-distance orbit is one, so every orbit-invariant selector is constant across representatives; raw-coordinate anchors distinguish representatives only after importing an external unit/gauge premise.  "
        "Therefore this exports no scale-invariant beta selector, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2654/S1604 scale-invariant beta selector Ltotal guard",
        "## P2654/S1604 scale-invariant beta selector Ltotal guard\n\n"
        "`P2654/S1604` does not re-enable `L_total`: invariant denominator functionals cannot select a unique beta representative, and raw-coordinate anchors are gauge/unit premises until sourced by typed nadsoliton metric/UV dynamics.  A variational damping coefficient still requires a sourced UV unit plus an independent beta-selecting operator/conservation identity, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
