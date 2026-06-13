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
OUT = GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json"
MD = GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.md"

ETA = 9.0 / 5.0
SAMPLE_BETAS = [0.01, 0.09270338861541028, 1.0, 1.1473957999384183, 4.0]
SAMPLE_DISTANCES = [1.0, 2.0, 7.0, 12.0]

SOURCES = {
    "P2647_HOLDOUT_HARNESS": GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json",
    "P2648_MARGIN_RULE": GEN / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.json",
    "P2649_BETA_ROUTE_MATRIX": GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json",
    "P2650_CANONICAL_UNIT_NO_GO": GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_length_uv_unit_exported",
    "target_independent_beta_identity_exported",
    "beta_source_exported",
    "legacy_role_transfer_revalidated",
    "empirical_holdout_promoted_to_source",
    "bridge_completion_exported",
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
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean",
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
        "gauge_fixed_beta_content": (
            "gauge-fixed|gauge fixed|working normalization|normalization representative|bare beta|"
            "beta=1.*working|normalization orbit|canonical length"
        ),
        "rescaling_invariant_content": (
            "d'=a\\*d|beta -> beta/a\\^eta|coordinate rescaling|invariant|beta\\*d\\^eta|"
            "tail ratio|denominator orbit"
        ),
        "role_boundary_content": (
            "role-transfer|legacy role|inverse-hierarchy|modified/compressed successor|"
            "compression/locality-bias|empirical holdout"
        ),
        "source_nonclosure_content": (
            "target-independent beta|beta source|source theorem|role-bearing L_total|QW-2191|"
            "ToE closure|bridge completion"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for beta=1 gauge-fixed working-normalization contract, not packet-name search",
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


def gauge_transform(beta: float, d: float, scale: float) -> tuple[float, float]:
    return beta / (scale**ETA), scale * d


def beta_one_scale(beta: float) -> float:
    return beta ** (1.0 / ETA)


def orbit_invariance_witness() -> dict[str, Any]:
    rows = []
    for beta in SAMPLE_BETAS:
        scale = beta_one_scale(beta)
        for d in SAMPLE_DISTANCES:
            beta_prime, d_prime = gauge_transform(beta, d, scale)
            original = envelope(beta, d)
            transformed = envelope(beta_prime, d_prime)
            rows.append({
                "original_beta": beta,
                "original_d": d,
                "scale_a_to_beta_prime_one": scale,
                "beta_prime": beta_prime,
                "d_prime": d_prime,
                "invariant_x_beta_d_eta": beta * d**ETA,
                "original_envelope": original,
                "transformed_envelope": transformed,
                "absolute_error": abs(original - transformed),
            })
    return {
        "group_action": "d' = a*d, beta' = beta/a^eta preserves x=beta*d^eta and the envelope 1/(1+x).",
        "beta_one_gauge": "For beta>0 choose a=beta^(1/eta), giving beta'=1 and d'=beta^(1/eta)*d.",
        "rows": rows,
        "max_envelope_error": max(row["absolute_error"] for row in rows),
        "all_rows_invariant_to_roundoff": all(row["absolute_error"] < 1e-15 for row in rows),
    }


def tail_ratio_gauge_witness() -> dict[str, Any]:
    pairs = [(1.0, 7.0), (1.0, 12.0), (2.0, 8.0), (3.0, 9.0)]
    rows = []
    for beta in SAMPLE_BETAS:
        scale = beta_one_scale(beta)
        for near, far in pairs:
            original_ratio = envelope(beta, far) / envelope(beta, near)
            beta_prime, near_prime = gauge_transform(beta, near, scale)
            _, far_prime = gauge_transform(beta, far, scale)
            transformed_ratio = envelope(beta_prime, far_prime) / envelope(beta_prime, near_prime)
            raw_strict_ratio_without_distance_transform = envelope(1.0, far) / envelope(1.0, near)
            rows.append({
                "original_beta": beta,
                "near": near,
                "far": far,
                "scale_a_to_beta_prime_one": scale,
                "near_prime": near_prime,
                "far_prime": far_prime,
                "ratio_original": original_ratio,
                "ratio_in_beta_one_gauge_with_scaled_distances": transformed_ratio,
                "ratio_if_beta_set_to_one_but_raw_distances_not_scaled": raw_strict_ratio_without_distance_transform,
                "gauge_respecting_error": abs(original_ratio - transformed_ratio),
                "raw_distance_substitution_error": abs(original_ratio - raw_strict_ratio_without_distance_transform),
            })
    return {
        "rule": "Tail ratios are gauge-respecting only when the distance coordinates are scaled with beta; setting beta=1 while leaving raw distances unchanged is a new gauge choice, not an invariant operation.",
        "rows": rows,
        "max_gauge_respecting_error": max(row["gauge_respecting_error"] for row in rows),
        "max_raw_distance_substitution_error": max(row["raw_distance_substitution_error"] for row in rows),
        "raw_substitution_has_nonzero_errors": any(row["raw_distance_substitution_error"] > 1e-9 for row in rows if abs(row["original_beta"] - 1.0) > 1e-12),
    }


def claim_boundary_matrix() -> list[dict[str, Any]]:
    rows = [
        {
            "claim": "strict_beta_equals_one_working_gauge",
            "allowed_under_contract": True,
            "status": "ALLOWED_AS_GAUGE_FIXED_WORKING_NORMALIZATION",
            "required_label": "must say beta=1 is a chosen representative of the damping orbit",
        },
        {
            "claim": "target_independent_beta_source",
            "allowed_under_contract": False,
            "status": "BLOCKED_BY_P2649_P2650",
            "required_label": "needs canonical length/UV unit plus source identity",
        },
        {
            "claim": "modified_compressed_successor_semantics",
            "allowed_under_contract": True,
            "status": "ALLOWED_AS_GAUGE_DECLARED_DESCRIPTION",
            "required_label": "must declare the beta=1 gauge and raw distance convention",
        },
        {
            "claim": "unchanged_legacy_inverse_hierarchy_transfer",
            "allowed_under_contract": False,
            "status": "REJECTED",
            "required_label": "P2643-P2645 rejection remains active",
        },
        {
            "claim": "blind_holdout_empirical_support",
            "allowed_under_contract": True,
            "status": "ALLOWED_ONLY_IF_DATA_GAUGE_AND_UNIT_MAP_ARE_DECLARED",
            "required_label": "cannot replace beta source theorem",
        },
        {
            "claim": "role_bearing_ltotal_damping_term",
            "allowed_under_contract": False,
            "status": "BLOCKED",
            "required_label": "needs beta source, typed unit, role-transfer, and selector/source discharge",
        },
    ]
    return rows


def upstream_consistency() -> dict[str, Any]:
    p2649 = load_json(SOURCES["P2649_BETA_ROUTE_MATRIX"])
    p2650 = load_json(SOURCES["P2650_CANONICAL_UNIT_NO_GO"])
    p2647 = load_json(SOURCES["P2647_HOLDOUT_HARNESS"])
    p2648 = load_json(SOURCES["P2648_MARGIN_RULE"])
    return {
        "p2649_no_beta_source_routes_pass": p2649.get("closure_decision", {}).get("passing_beta_source_routes") == [],
        "p2650_no_canonical_unit_candidates_pass": p2650.get("closure_decision", {}).get("passing_candidates") == [],
        "p2650_recommends_gauge_fixed_working_normalization": "gauge-fixed working normalization" in json.dumps(p2650, sort_keys=True),
        "p2647_harness_ready_but_no_real_holdout": p2647.get("closure_decision", {}).get("real_blind_holdout_executed") is False,
        "p2648_margin_rule_ready_but_no_blind_data": p2648.get("closure_decision", {}).get("real_blind_holdout_executed") is False,
    }


def closure_decision(boundary: list[dict[str, Any]]) -> dict[str, Any]:
    blocked = [row["claim"] for row in boundary if not row["allowed_under_contract"]]
    allowed = [row["claim"] for row in boundary if row["allowed_under_contract"]]
    return {
        "decision": "BETA_ONE_GAUGE_FIXED_WORKING_NORMALIZATION_CONTRACT_ONLY__NO_SOURCE_NO_LTOTAL_NO_TOE",
        "allowed_claims": allowed,
        "blocked_claims": blocked,
        "professorial_verdict": (
            "P2651 is the honest fallback after P2649-P2650: beta=1 may be used as a declared gauge-fixed working normalization, "
            "because every positive beta lies on a rescaling orbit with a beta=1 representative.  The contract prevents that representative from being promoted to a beta source, legacy role transfer, or role-bearing dynamics."
        ),
        "professorial_closure_path": [
            "Use beta=1 only with an explicit gauge declaration and transformed distance convention.",
            "For empirical P2647/P2648 runs, require a data-unit map into the declared beta=1 gauge before comparing tail ratios.",
            "For theory closure, the missing next atom is still a typed nadsoliton metric/UV source, not a larger numerical scan.",
            "Rerun role-transfer only after beta source, bridge completion, and selector/source blockers are independently discharged.",
        ],
        "next_honest_step": (
            "Adopt the P2651 gauge-fixed contract in downstream packets, then either build a typed nadsoliton metric/UV source theorem or execute P2647/P2648 with a declared data-unit map as empirical support only."
        ),
        "beta_one_working_gauge_allowed": True,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    lines = [
        "# P2651/S1601 beta=1 gauge-fixed working-normalization contract",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps gauge-fixed beta, rescaling invariants, role boundaries, and source nonclosure before adding the contract.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Orbit invariance",
        "",
        payload["orbit_invariance_witness"]["group_action"],
        payload["orbit_invariance_witness"]["beta_one_gauge"],
        f"Max envelope error on audited grid: `{payload['orbit_invariance_witness']['max_envelope_error']:.3e}`.",
        "",
        "## Tail-ratio gauge warning",
        "",
        payload["tail_ratio_gauge_witness"]["rule"],
        f"Max gauge-respecting ratio error: `{payload['tail_ratio_gauge_witness']['max_gauge_respecting_error']:.3e}`.",
        f"Max raw-distance substitution error: `{payload['tail_ratio_gauge_witness']['max_raw_distance_substitution_error']:.12f}`.",
        "",
        "## Claim boundary matrix",
        "",
        "| claim | allowed? | status | required label |",
        "| --- | ---: | --- | --- |",
    ])
    for row in payload["claim_boundary_matrix"]:
        lines.append(f"| `{row['claim']}` | `{row['allowed_under_contract']}` | `{row['status']}` | {row['required_label']} |")
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Beta=1 working gauge allowed? `{decision['beta_one_working_gauge_allowed']}`.",
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
    boundary = claim_boundary_matrix()
    payload: dict[str, Any] = {
        "status": "P2651_BETA_ONE_GAUGE_FIXED_WORKING_NORMALIZATION_CONTRACT_NO_SOURCE_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCES.items()},
        "upstream_consistency": upstream_consistency(),
        "orbit_invariance_witness": orbit_invariance_witness(),
        "tail_ratio_gauge_witness": tail_ratio_gauge_witness(),
        "claim_boundary_matrix": boundary,
        "closure_decision": closure_decision(boundary),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCES["STRICT_EQUATION_SHEET"],
        "P2651/S1601 beta=1 gauge-fixed working-normalization guard",
        "## P2651/S1601 beta=1 gauge-fixed working-normalization guard\n\n"
        "`P2651/S1601` formalizes the honest fallback after P2649-P2650: `beta=1` is allowed as a declared gauge-fixed working normalization because `d'=a*d`, `beta'=beta/a^eta` preserves `beta*d^eta` and every positive beta has a `beta'=1` representative.  Tail ratios remain gauge-respecting only when distances are transformed with the unit map; setting `beta=1` while leaving raw distances unchanged is a new gauge choice, not an invariant operation.  This contract allows modified/compressed successor semantics and empirical holdout support only with explicit gauge/unit declarations, but exports no beta source, legacy role transfer, bridge completion, `QW-2191` discharge, role-bearing `L_total`, or ToE closure.\n",
    )
    append_once(
        SOURCES["STRICT_LAGRANGIAN_DRAFT"],
        "P2651/S1601 beta=1 gauge-fixed Ltotal guard",
        "## P2651/S1601 beta=1 gauge-fixed Ltotal guard\n\n"
        "`P2651/S1601` does not re-enable `L_total`: it permits `beta=1` only as a gauge-fixed working normalization with explicit unit-map bookkeeping, not as a sourced variational damping coefficient.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
