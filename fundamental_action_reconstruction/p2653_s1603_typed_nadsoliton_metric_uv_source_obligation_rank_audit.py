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
OUT = GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json"
MD = GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.md"

P2649 = GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json"
P2650 = GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json"
P2651 = GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json"
P2652 = GEN / "p2652_s1602_beta_one_gauge_unit_map_validator_and_holdout_covariance_firewall.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

ETA = 9.0 / 5.0
AUDITED_BETAS = [0.01, 0.09270338861541028, 1.0, 1.1473957999384183, 4.0]
AUDITED_DISTANCES = [1.0, 2.0, 7.0, 12.0]
REQUIRED_TYPED_UNIT_ATOMS = [
    "typed_nadsoliton_state_space",
    "positive_metric_distance_functional",
    "uv_unit_selected_by_nadsoliton_dynamics",
    "unit_selection_independent_of_strict_kernel_target",
    "unit_selection_independent_of_empirical_holdout_fit",
    "scale_orbit_quotient_discharge",
    "dimensionless_conservation_or_operator_identity",
    "unique_beta_one_solution_after_unit_fixing",
]

NEGATIVE_EXPORT_FLAGS = [
    "typed_metric_uv_source_theorem_exported",
    "canonical_unit_exported",
    "scale_orbit_discharged",
    "target_independent_beta_source_exported",
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
        "typed_metric_uv_content": (
            "typed nadsoliton metric|metric/UV|UV unit|canonical length|distance functional|"
            "nadsoliton dynamics.*unit|state space"
        ),
        "scale_orbit_obstruction_content": (
            "scale orbit|normalization orbit|gauge-fixed|gauge fixed|beta -> beta/a\\^eta|"
            "unit map|raw-to-beta|scale convention"
        ),
        "operator_identity_content": (
            "operator identity|conservation identity|dimensionless conservation|source theorem|"
            "unique beta|target-independent beta"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|"
            "empirical confirmation"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for typed nadsoliton metric/UV source obligations, not packet-name search",
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


def scale_to_target_beta(beta: float, target_beta: float) -> float:
    return (beta / target_beta) ** (1.0 / ETA)


def scale_orbit_indeterminacy_witness() -> dict[str, Any]:
    rows = []
    for beta in AUDITED_BETAS:
        for target_beta in AUDITED_BETAS:
            scale = scale_to_target_beta(beta, target_beta)
            max_error = 0.0
            for d in AUDITED_DISTANCES:
                original = envelope(beta, d)
                transformed = envelope(target_beta, scale * d)
                max_error = max(max_error, abs(original - transformed))
            rows.append({
                "source_beta": beta,
                "target_beta_representative": target_beta,
                "scale_a": scale,
                "max_envelope_error_on_grid": max_error,
            })
    return {
        "statement": "For every positive beta and target representative beta*, d'=(beta/beta*)^(1/eta)d preserves beta*d^eta, so numeric beta is not fixed until a typed unit fixes the scale orbit.",
        "rows": rows,
        "max_error": max(row["max_envelope_error_on_grid"] for row in rows),
        "all_representatives_equivalent_to_roundoff": all(row["max_envelope_error_on_grid"] < 1e-15 for row in rows),
        "orbit_size_audited": len(rows),
    }


def route_atom_matrix() -> list[dict[str, Any]]:
    routes = [
        {
            "route": "dimensionless_coordinate_d_equals_1",
            "available_atoms": {"typed_nadsoliton_state_space"},
            "failure_mode": "coordinate convention supplies no nadsoliton-selected UV unit and leaves scale orbit free",
        },
        {
            "route": "legacy_beta_tors_or_alpha_geo_unit",
            "available_atoms": {"typed_nadsoliton_state_space"},
            "failure_mode": "legacy numeric constants need bridge completion and role-transfer before they can type a strict metric unit",
        },
        {
            "route": "micro_zbeta_or_renormalization_unit",
            "available_atoms": {"typed_nadsoliton_state_space", "positive_metric_distance_functional"},
            "failure_mode": "micro/strict mismatch and normalization gauge are not discharged by a unique UV selector",
        },
        {
            "route": "empirical_tail_ratio_unit_calibration",
            "available_atoms": {"positive_metric_distance_functional", "unit_selection_independent_of_strict_kernel_target"},
            "failure_mode": "empirical calibration can support compression but is still a holdout/target route, not a pre-kernel unit source",
        },
        {
            "route": "p2652_precommitted_unit_map_validator",
            "available_atoms": {
                "typed_nadsoliton_state_space",
                "positive_metric_distance_functional",
                "unit_selection_independent_of_empirical_holdout_fit",
            },
            "failure_mode": "validator checks bookkeeping but does not export the source theorem that selects the unit map",
        },
        {
            "route": "hypothetical_valid_metric_uv_theorem",
            "available_atoms": set(REQUIRED_TYPED_UNIT_ATOMS),
            "failure_mode": "not currently in repo; included only as a theorem-target specification",
        },
    ]
    matrix = []
    for route in routes:
        available = route["available_atoms"]
        missing = [atom for atom in REQUIRED_TYPED_UNIT_ATOMS if atom not in available]
        matrix.append({
            "route": route["route"],
            "available_atoms": sorted(available),
            "missing_atoms": missing,
            "available_atom_count": len(available),
            "required_atom_count": len(REQUIRED_TYPED_UNIT_ATOMS),
            "passes_typed_metric_uv_source": not missing,
            "failure_mode": route["failure_mode"],
        })
    return matrix


def obligation_rank_audit(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    current_rows = [row for row in matrix if row["route"] != "hypothetical_valid_metric_uv_theorem"]
    covered = sorted(set().union(*(set(row["available_atoms"]) for row in current_rows))) if current_rows else []
    missing_globally = [atom for atom in REQUIRED_TYPED_UNIT_ATOMS if atom not in covered]
    strongest_current = max(current_rows, key=lambda row: row["available_atom_count"])
    return {
        "required_atoms": REQUIRED_TYPED_UNIT_ATOMS,
        "currently_covered_atoms_union": covered,
        "currently_missing_atoms_union": missing_globally,
        "strongest_current_route": strongest_current["route"],
        "strongest_current_route_missing_atoms": strongest_current["missing_atoms"],
        "current_routes_passing": [row["route"] for row in current_rows if row["passes_typed_metric_uv_source"]],
        "hypothetical_theorem_target_passing": [row["route"] for row in matrix if row["passes_typed_metric_uv_source"]],
        "current_atom_coverage_fraction": len(covered) / len(REQUIRED_TYPED_UNIT_ATOMS),
        "rank_defect_interpretation": "Current artifacts cover bookkeeping and some metric-like language, but miss the source atoms that choose a UV unit and a dimensionless operator/conservation identity before the strict target is read.",
    }


def upstream_consistency() -> dict[str, Any]:
    p2649 = load_json(P2649)
    p2650 = load_json(P2650)
    p2651 = load_json(P2651)
    p2652 = load_json(P2652)
    return {
        "p2649_beta_routes_blocked": p2649.get("closure_decision", {}).get("passing_beta_source_routes") == [],
        "p2650_no_canonical_unit_candidates_pass": p2650.get("closure_decision", {}).get("passing_candidates") == [],
        "p2651_beta_one_only_working_gauge": p2651.get("closure_decision", {}).get("beta_one_working_gauge_allowed") is True,
        "p2652_unit_map_validator_ready_not_source": p2652.get("closure_decision", {}).get("unit_map_source_theorem_exported_now") is False,
        "p2652_real_payload_not_loaded": p2652.get("closure_decision", {}).get("real_blind_holdout_payload_loaded") is False,
    }


def closure_decision(rank: dict[str, Any], orbit: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "TYPED_METRIC_UV_SOURCE_OBLIGATION_SPECIFIED_BUT_NOT_PROVED__SCALE_ORBIT_REMAINS_OPEN",
        "professorial_verdict": (
            "P2653 is the proof-facing audit after P2652: it states the typed nadsoliton metric/UV source theorem as explicit obligations and checks the current routes against them.  "
            "The scale-orbit calculation shows that any positive beta can be represented by any other positive beta after a distance rescaling, so beta=1 cannot become sourced without a prior unit selector.  "
            "Current artifacts validate bookkeeping, gauges, and empirical harness readiness, but they do not supply the typed UV unit or the independent operator/conservation identity that would uniquely select beta=1."
        ),
        "professorial_closure_path": [
            "Do not add larger numeric expression scans or more target-fitted unit maps; they cannot discharge the scale orbit.",
            "The next proof attempt must construct a typed nadsoliton state space, metric distance, and UV unit chosen by nadsoliton dynamics before reading K_strict_gate.",
            "Only after that unit is fixed should P2649 be rerun to test whether a dimensionless conservation/operator identity has beta=1 as a unique solution.",
            "Empirical P2652/P2647/P2648 remains useful only as compression support with a precommitted unit map, not as a beta-source theorem.",
        ],
        "next_honest_step": (
            "Try to prove the typed nadsoliton metric/UV source theorem by supplying the missing atoms `uv_unit_selected_by_nadsoliton_dynamics`, `scale_orbit_quotient_discharge`, and `dimensionless_conservation_or_operator_identity`; otherwise keep beta=1 as gauge-fixed working normalization and run real blinded payloads only through P2652/P2647/P2648 as support."
        ),
        "scale_orbit_equivalence_verified": orbit["all_representatives_equivalent_to_roundoff"],
        "current_routes_passing_typed_metric_uv_source": rank["current_routes_passing"],
        "typed_metric_uv_source_theorem_exported_now": False,
        "canonical_unit_exported_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    rank = payload["obligation_rank_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2653/S1603 typed nadsoliton metric/UV source obligation rank audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps typed metric/UV, scale-orbit obstruction, operator identity, and nonclosure guard content before adding the audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Scale-orbit witness",
        "",
        payload["scale_orbit_indeterminacy_witness"]["statement"],
        f"Audited beta-representative pairs: `{payload['scale_orbit_indeterminacy_witness']['orbit_size_audited']}`.",
        f"Max envelope error on audited grid: `{payload['scale_orbit_indeterminacy_witness']['max_error']:.3e}`.",
        "",
        "## Obligation matrix",
        "",
        "| route | passes? | available atoms | missing atoms |",
        "| --- | ---: | ---: | --- |",
    ])
    for row in payload["route_atom_matrix"]:
        lines.append(
            f"| `{row['route']}` | `{row['passes_typed_metric_uv_source']}` | "
            f"`{row['available_atom_count']}/{row['required_atom_count']}` | `{', '.join(row['missing_atoms'])}` |"
        )
    lines.extend([
        "",
        "## Rank audit",
        "",
        f"Current routes passing: `{rank['current_routes_passing']}`.",
        f"Current atom coverage fraction: `{rank['current_atom_coverage_fraction']:.3f}`.",
        rank["rank_defect_interpretation"],
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Typed metric/UV source theorem exported now? `{decision['typed_metric_uv_source_theorem_exported_now']}`.",
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
    matrix = route_atom_matrix()
    orbit = scale_orbit_indeterminacy_witness()
    rank = obligation_rank_audit(matrix)
    decision = closure_decision(rank, orbit)
    payload: dict[str, Any] = {
        "status": "P2653_TYPED_METRIC_UV_SOURCE_OBLIGATION_SPECIFIED_NOT_PROVED_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2649_BETA_ROUTE_MATRIX": sha256_file(P2649),
            "P2650_CANONICAL_UNIT_NO_GO": sha256_file(P2650),
            "P2651_GAUGE_CONTRACT": sha256_file(P2651),
            "P2652_UNIT_MAP_VALIDATOR": sha256_file(P2652),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "required_typed_metric_uv_atoms": REQUIRED_TYPED_UNIT_ATOMS,
        "scale_orbit_indeterminacy_witness": orbit,
        "route_atom_matrix": matrix,
        "obligation_rank_audit": rank,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2653/S1603 typed nadsoliton metric/UV source obligation guard",
        "## P2653/S1603 typed nadsoliton metric/UV source obligation guard\n\n"
        "`P2653/S1603` specifies the proof obligations for promoting `beta=1` beyond a gauge-fixed working normalization: a typed nadsoliton state space, metric distance, nadsoliton-selected UV unit, target/holdout independence, scale-orbit quotient discharge, and a dimensionless operator/conservation identity with unique `beta=1`.  "
        "The audited current routes cover only bookkeeping and partial metric language; the scale-orbit witness shows every positive beta remains representable by distance rescaling.  "
        "Thus this exports no typed metric/UV source theorem, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2653/S1603 typed metric UV Ltotal guard",
        "## P2653/S1603 typed metric UV Ltotal guard\n\n"
        "`P2653/S1603` does not re-enable `L_total`: it turns the missing typed nadsoliton metric/UV source into an explicit obligation matrix and verifies the scale-orbit obstruction remains open.  A variational damping coefficient still requires a sourced UV unit plus an independent beta-selecting operator/conservation identity, real blinded support, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
