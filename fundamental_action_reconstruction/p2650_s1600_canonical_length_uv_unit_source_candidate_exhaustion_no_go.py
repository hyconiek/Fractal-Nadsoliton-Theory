#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json"
MD = GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.md"

ETA = 9.0 / 5.0
ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
BETA_CRIT = 0.09270338861541028
MICRO_BETA_MEDIAN = 1.1473957999384183
STRICT_BETA = 1.0

SOURCES = {
    "P2629_ZBETA_GAUGE": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2631_BETA_CRITICALITY": GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json",
    "P2643_BETA_THRESHOLD": GEN / "p2643_s1593_inverse_hierarchy_beta_threshold_role_rejection_certificate.json",
    "P2649_BETA_ROUTE_MATRIX": GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_length_uv_unit_exported",
    "target_independent_beta_identity_exported",
    "normalization_gauge_fixed",
    "micro_strict_mismatch_removed",
    "legacy_beta_tors_promoted_to_strict_source",
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
        "canonical_length_uv_content": (
            "canonical length|canonical UV|UV unit|length/UV|normalization gauge|scale convention|"
            "bare beta|beta=1.*normalization|canonical unit"
        ),
        "nadsoliton_source_content": (
            "nadsoliton dynamics|source theorem|target-independent|operator identity|conservation constant|"
            "dimensionless conservation|first principles"
        ),
        "legacy_strict_unit_bridge_content": (
            "beta_tors|alpha_geo|legacy.*strict|completion bridge|bridge completion|"
            "D_f|eta=9/5|strict beta"
        ),
        "micro_mismatch_content": (
            "micro beta|Z_beta|beta_micro|beta_strict|normalization-invariant|micro/strict mismatch|"
            "renormalization constants"
        ),
        "closure_guard_content": (
            "role-transfer|role-bearing L_total|QW-2191|ToE closure|empirical confirmation|"
            "positive_beta_renormalization_source"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for canonical length/UV unit source candidates, not packet-name search",
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


def canonical_unit_candidates() -> list[dict[str, Any]]:
    base_candidates = [
        {
            "candidate": "dimensionless_domain_unit_d_equals_1",
            "proposed_unit_a": 1.0,
            "source_claim": "take the audited strict coordinate unit d=1 as canonical",
            "required_atoms": {"nadsoliton_metric_unit_source", "uv_normalization_gauge_fixed", "not_chosen_from_strict_target"},
            "available_atoms": {"coordinate_unit_available"},
        },
        {
            "candidate": "legacy_beta_tors_unit",
            "proposed_unit_a": BETA_TORS ** (1.0 / ETA),
            "source_claim": "use legacy beta_tors to define the strict length unit",
            "required_atoms": {"beta_tors_to_strict_beta_source_map", "role_transfer_theorem", "uv_normalization_gauge_fixed"},
            "available_atoms": {"legacy_beta_tors_numeric_available"},
        },
        {
            "candidate": "alpha_geo_unit",
            "proposed_unit_a": ALPHA_GEO ** (1.0 / ETA),
            "source_claim": "use alpha_geo information-geometry amplitude as length/UV unit",
            "required_atoms": {"alpha_geo_to_damping_operator_source", "dimensionally_typed_length_map", "uv_normalization_gauge_fixed"},
            "available_atoms": {"alpha_geo_numeric_available"},
        },
        {
            "candidate": "df_eta_unit",
            "proposed_unit_a": ETA ** (1.0 / ETA),
            "source_claim": "use D_f=eta=9/5 as canonical fractal length unit",
            "required_atoms": {"fractal_exponent_to_length_unit_theorem", "damping_coefficient_source", "uv_normalization_gauge_fixed"},
            "available_atoms": {"df_equals_eta_available"},
        },
        {
            "candidate": "inverse_hierarchy_threshold_unit",
            "proposed_unit_a": BETA_CRIT ** (1.0 / ETA),
            "source_claim": "use the P2643 inverse-hierarchy beta threshold as canonical unit",
            "required_atoms": {"threshold_promoted_to_source_theorem", "unchanged_inverse_hierarchy_role_revalidated", "uv_normalization_gauge_fixed"},
            "available_atoms": {"threshold_numeric_available"},
        },
        {
            "candidate": "micro_beta_median_unit",
            "proposed_unit_a": MICRO_BETA_MEDIAN ** (1.0 / ETA),
            "source_claim": "use the target-blind micro beta median to define the strict unit",
            "required_atoms": {"micro_strict_mismatch_removed", "positive_zbeta_source_exported", "uv_normalization_gauge_fixed"},
            "available_atoms": {"micro_beta_median_numeric_available"},
        },
        {
            "candidate": "empirical_tail_ratio_unit",
            "proposed_unit_a": 1.0,
            "source_claim": "let future blind tail-ratio data set the strict beta unit",
            "required_atoms": {"real_blind_holdout_passed", "empirical_holdout_can_replace_source_theorem", "uv_normalization_gauge_fixed"},
            "available_atoms": {"p2647_p2648_protocol_available"},
        },
    ]
    rows = []
    for row in base_candidates:
        atom_truth = {atom: atom in row["available_atoms"] for atom in sorted(row["available_atoms"] | row["required_atoms"])}
        missing = sorted(atom for atom in row["required_atoms"] if not atom_truth.get(atom, False))
        beta_prime_from_strict_one = STRICT_BETA / (row["proposed_unit_a"] ** ETA)
        rows.append({
            "candidate": row["candidate"],
            "source_claim": row["source_claim"],
            "proposed_unit_a": row["proposed_unit_a"],
            "beta_prime_if_strict_beta_one_is_rescaled_by_a": beta_prime_from_strict_one,
            "required_atoms": sorted(row["required_atoms"]),
            "available_atoms": sorted(row["available_atoms"]),
            "missing_required_atoms": missing,
            "passes_as_canonical_length_uv_source": not missing,
            "rejection_reason": "PASS" if not missing else "missing " + ", ".join(missing),
        })
    return rows


def finite_expression_scan() -> dict[str, Any]:
    constants = {
        "alpha_geo": ALPHA_GEO,
        "D_f_eta": ETA,
        "beta_tors": BETA_TORS,
        "beta_crit": BETA_CRIT,
        "micro_beta_median": MICRO_BETA_MEDIAN,
    }
    exponents = [-1.0, -1.0 / ETA, -0.5, 0.5, 1.0 / ETA, 1.0]
    rows = []
    for (name_a, value_a), (name_b, value_b), exp_a, exp_b in itertools.product(constants.items(), constants.items(), exponents, exponents):
        unit = (value_a ** exp_a) * (value_b ** exp_b)
        beta_prime = STRICT_BETA / (unit ** ETA)
        near_one = abs(beta_prime - 1.0)
        if near_one < 0.05:
            rows.append({
                "expression": f"{name_a}^{exp_a:.6g} * {name_b}^{exp_b:.6g}",
                "unit_a": unit,
                "beta_prime_from_strict_one": beta_prime,
                "distance_to_beta_one": near_one,
                "source_status": "numeric expression only; no typed nadsoliton length-source theorem",
            })
    rows = sorted(rows, key=lambda row: row["distance_to_beta_one"])[:20]
    return {
        "grammar": "two-factor products of audited constants with exponents {-1, -1/eta, -1/2, 1/2, 1/eta, 1}",
        "near_beta_one_rows": rows,
        "near_rows_count": len(rows),
        "minimum_distance_to_beta_one": rows[0]["distance_to_beta_one"] if rows else None,
        "exports_canonical_unit": False,
        "verdict": "Finite numerology can find units close to beta'=1, but every expression is still untyped and source-free unless a nadsoliton length/UV theorem selects it before comparison to the strict target.",
    }


def upstream_consistency() -> dict[str, Any]:
    p2629 = load_json(SOURCES["P2629_ZBETA_GAUGE"])
    p2643 = load_json(SOURCES["P2643_BETA_THRESHOLD"])
    p2649 = load_json(SOURCES["P2649_BETA_ROUTE_MATRIX"])
    return {
        "p2629_normalization_gauge_unfixed": not bool(p2629.get("exact_source_gate", {}).get("gates", {}).get("normalization_invariant_mismatch_removed", False)),
        "p2629_micro_over_strict_ratio": p2629.get("normalization_orbit_certificate", {}).get("invariant_ratio"),
        "p2643_beta_crit_available": p2643.get("beta_threshold_theorem", {}).get("beta_critical_exact_role_boundary"),
        "p2649_no_beta_source_routes_pass": p2649.get("closure_decision", {}).get("passing_beta_source_routes") == [],
        "p2649_demands_canonical_length_or_conservation_identity": "canonical-length" in json.dumps(p2649, sort_keys=True) or "canonical length" in json.dumps(p2649, sort_keys=True),
    }


def closure_decision(candidate_rows: list[dict[str, Any]], expression_scan: dict[str, Any]) -> dict[str, Any]:
    passing = [row["candidate"] for row in candidate_rows if row["passes_as_canonical_length_uv_source"]]
    return {
        "decision": "NO_CANONICAL_LENGTH_UV_UNIT_SOURCE_CANDIDATE_PASSES__BETA_SOURCE_REMAINS_BLOCKED",
        "passing_candidates": passing,
        "professorial_verdict": (
            "P2650 audits the exact atom P2649 said was needed: a canonical length/UV unit independent of the strict target.  "
            "The audited candidates either choose a coordinate convention, reuse legacy beta_tors without a role-transfer theorem, reuse micro/empirical targets, or form untyped numerology.  "
            "Therefore beta=1 remains a robust working normalization/compression parameter, not sourced strict dynamics."
        ),
        "professorial_closure_path": [
            "Do not search larger numeric expression grammars for beta=1; they only create more unsourced unit choices.",
            "A valid next theorem must construct a typed nadsoliton metric/UV unit before reading the strict denominator.",
            "After such a unit exists, rerun P2649: only then can a dimensionless conservation/operator identity be tested for beta=1 without circularity.",
            "If no typed unit theorem emerges, keep the empirical P2647/P2648 path as support only and label beta sourcehood blocked.",
        ],
        "next_honest_step": (
            "Build a typed nadsoliton metric/UV-unit theorem that fixes the coordinate d=1 independently of K_strict_gate; if that cannot be done, formalize beta=1 as a gauge-fixed working normalization and keep role-bearing L_total closed."
        ),
        "finite_expression_scan_exports_unit": expression_scan["exports_canonical_unit"],
        "canonical_length_uv_unit_exported_now": False,
        "beta_source_exported_now": False,
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    lines = [
        "# P2650/S1600 canonical length/UV unit source candidate exhaustion no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps canonical length/UV, nadsoliton source, legacy-strict unit bridge, micro mismatch, and closure guard content before adding the audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Candidate matrix",
        "",
        "| candidate | passes? | proposed unit a | beta' from strict beta=1 | missing atoms |",
        "| --- | ---: | ---: | ---: | --- |",
    ])
    for row in payload["canonical_unit_candidates"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['passes_as_canonical_length_uv_source']}` | "
            f"`{row['proposed_unit_a']:.12g}` | `{row['beta_prime_if_strict_beta_one_is_rescaled_by_a']:.12g}` | "
            f"`{', '.join(row['missing_required_atoms'])}` |"
        )
    lines.extend([
        "",
        "## Finite expression scan",
        "",
        payload["finite_expression_scan"]["grammar"],
        payload["finite_expression_scan"]["verdict"],
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing candidates: `{decision['passing_candidates']}`.",
        f"Canonical length/UV exported now? `{decision['canonical_length_uv_unit_exported_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
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
    candidates = canonical_unit_candidates()
    scan = finite_expression_scan()
    decision = closure_decision(candidates, scan)
    payload: dict[str, Any] = {
        "status": "P2650_CANONICAL_LENGTH_UV_UNIT_SOURCE_EXHAUSTION_NO_GO_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCES.items()},
        "upstream_consistency": upstream_consistency(),
        "canonical_unit_candidates": candidates,
        "finite_expression_scan": scan,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCES["STRICT_EQUATION_SHEET"],
        "P2650/S1600 canonical length/UV unit source guard",
        "## P2650/S1600 canonical length/UV unit source guard\n\n"
        "`P2650/S1600` audits the canonical length/UV unit atom required by P2649.  The audited candidates `d=1`, `beta_tors`, `alpha_geo`, `D_f=eta`, the inverse-hierarchy threshold, the micro beta median, and future empirical tail-ratio calibration all fail as target-independent unit sources because they lack typed nadsoliton metric/UV source atoms, role-transfer theorems, or normalization-gauge discharge.  A finite expression scan over audited constants can generate near-unit numerology but still exports no typed unit theorem.  Thus `beta=1` remains a robust working normalization/compression parameter, not a sourced strict beta; bridge completion, role transfer, `QW-2191`, role-bearing `L_total`, and ToE closure remain closed.\n",
    )
    append_once(
        SOURCES["STRICT_LAGRANGIAN_DRAFT"],
        "P2650/S1600 canonical length/UV unit Ltotal guard",
        "## P2650/S1600 canonical length/UV unit Ltotal guard\n\n"
        "`P2650/S1600` does not re-enable `L_total`: no audited canonical length/UV candidate supplies the typed nadsoliton unit needed to turn `beta=1` from gauge-fixed working normalization into a variational source parameter.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
