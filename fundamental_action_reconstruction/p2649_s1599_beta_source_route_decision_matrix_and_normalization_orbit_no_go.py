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
OUT = GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json"
MD = GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.md"

ETA = 9.0 / 5.0
AUDITED_BETAS = [0.01, 0.0927, 1.0, 1.1473958, 4.0]
TAIL_PAIRS = [(1, 7), (1, 12), (2, 8), (3, 9)]

SOURCES = {
    "P2629_ZBETA_GAUGE": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2630_BETA_SOURCE_SEPARATION": GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json",
    "P2631_BETA_CRITICALITY": GEN / "p2631_s1581_neural_information_flux_beta_criticality_audit.json",
    "P2645_ROLE_MATRIX": GEN / "p2645_s1595_role_transfer_matrix_and_closure_route_rerun.json",
    "P2647_HOLDOUT_HARNESS": GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json",
    "P2648_MARGIN_RULE": GEN / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.json",
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "STRICT_LAGRANGIAN_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
    "target_independent_beta_identity_exported",
    "normalization_gauge_fixed",
    "positive_zbeta_source_exported",
    "micro_strict_mismatch_removed",
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
        "beta_source_identity_content": (
            "target-independent beta|beta-source|beta source|beta=1.*first principles|"
            "positive_beta_renormalization_source|strict beta source|beta identity"
        ),
        "normalization_orbit_content": (
            "normalization gauge|UV normalization|length/UV|scale convention|d -> a\\*d|"
            "beta -> beta/a\\^eta|bare beta|canonical length"
        ),
        "micro_zbeta_mismatch_content": (
            "Z_beta|beta_micro|beta_strict|micro/strict mismatch|beta_micro/beta_strict|"
            "target-blind micro|renormalization constants"
        ),
        "empirical_not_source_content": (
            "blind holdout|frozen-kernel|statistical margin|empirical confirmation|"
            "synthetic fixture|standard error|Bonferroni"
        ),
        "closure_guard_content": (
            "role-transfer|bridge completion|role-bearing L_total|QW-2191|ToE closure|"
            "source theorem"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for beta-source route decision matrix and normalization-orbit no-go, not packet-name search",
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


def tail_ratio(beta: float, near: float, far: float) -> float:
    return (1.0 + beta * near**ETA) / (1.0 + beta * far**ETA)


def beta_from_tail_ratio(q: float, near: float, far: float) -> float | None:
    denominator = q * far**ETA - near**ETA
    if denominator <= 0:
        return None
    beta = (1.0 - q) / denominator
    return beta if beta > 0 else None


def normalization_orbit_witness() -> dict[str, Any]:
    rows = []
    for beta in AUDITED_BETAS:
        scale_to_beta_one = beta ** (1.0 / ETA)
        beta_after_scale = beta / (scale_to_beta_one**ETA)
        invariant = beta / (scale_to_beta_one**ETA)
        rows.append({
            "original_beta": beta,
            "length_rescaling_a_that_sets_beta_prime_to_1": scale_to_beta_one,
            "beta_prime_after_rescaling": beta_after_scale,
            "denominator_invariant_check_beta_over_a_eta": invariant,
            "bare_beta_one_reached_by_rescaling": abs(beta_after_scale - 1.0) < 1e-12,
        })
    return {
        "group_action": "For the envelope denominator 1+beta*d^eta, coordinate rescaling d'=a*d sends beta to beta'=beta/a^eta.",
        "theorem": "For every beta>0 there exists a positive a=beta^(1/eta) that rewrites the same denominator orbit with beta'=1.",
        "consequence": "The bare numeral beta=1 is an orbit representative, not a source theorem, until a canonical length/UV unit is independently fixed.",
        "rows": rows,
        "all_positive_betas_gauge_to_one": all(row["bare_beta_one_reached_by_rescaling"] for row in rows),
        "exports_beta_source": False,
    }


def tail_ratio_inversion_witness() -> dict[str, Any]:
    rows = []
    for near, far in TAIL_PAIRS:
        q_strict = tail_ratio(1.0, near, far)
        beta_recovered = beta_from_tail_ratio(q_strict, near, far)
        q_legacy = tail_ratio(0.01, near, far)
        beta_legacy_recovered = beta_from_tail_ratio(q_legacy, near, far)
        rows.append({
            "near": near,
            "far": far,
            "strict_ratio_q_beta_1": q_strict,
            "beta_recovered_from_strict_ratio": beta_recovered,
            "legacy_ratio_q_beta_0_01": q_legacy,
            "beta_recovered_from_legacy_ratio": beta_legacy_recovered,
            "lesson": "a chosen tail ratio fixes beta algebraically; it is a target datum unless independently sourced",
        })
    return {
        "formula": "For R(a,b)=q=(1+beta*a^eta)/(1+beta*b^eta), beta=(1-q)/(q*b^eta-a^eta) when the denominator is positive.",
        "consequence": "A single observed or imposed tail ratio can recover beta=1, but only by encoding beta=1 in q; this is an empirical/calibration target, not target-independent sourcehood.",
        "rows": rows,
        "all_strict_ratios_recover_beta_one": all(abs(row["beta_recovered_from_strict_ratio"] - 1.0) < 1e-12 for row in rows if row["beta_recovered_from_strict_ratio"] is not None),
        "exports_beta_source": False,
    }


def upstream_atoms() -> dict[str, bool]:
    p2629 = load_json(SOURCES["P2629_ZBETA_GAUGE"])
    p2630 = load_json(SOURCES["P2630_BETA_SOURCE_SEPARATION"])
    p2631 = load_json(SOURCES["P2631_BETA_CRITICALITY"])
    p2645 = load_json(SOURCES["P2645_ROLE_MATRIX"])
    p2647 = load_json(SOURCES["P2647_HOLDOUT_HARNESS"])
    p2648 = load_json(SOURCES["P2648_MARGIN_RULE"])
    return {
        "normalization_gauge_fixed": False,
        "canonical_length_unit_exported": False,
        "target_independent_conservation_constant_exported": False,
        "unique_beta_one_critical_extremum_exported": False,
        "positive_zbeta_source_exported": bool(p2629.get("exact_source_gate", {}).get("accepts_positive_beta_renormalization_source", False)),
        "bridge_zbeta_source_accepted": bool(p2630.get("bridge_source_truth_table", {}).get("current_accepts_bridge_positive_zbeta_source", False)),
        "p2631_exports_beta_identity": bool(p2631.get("source_acceptance", {}).get("accepts_information_conservation_beta_identity", False)),
        "role_matrix_beta_source_passes": "strict_beta_source_role" in json.dumps(p2645, sort_keys=True) and "target_independent_beta_identity" not in json.dumps(p2645.get("role_transfer_matrix", []), sort_keys=True),
        "holdout_harness_ready": bool(p2647.get("fake_pass_firewall", {}).get("firewall_passes", False)),
        "holdout_real_blind_data_executed": bool(p2647.get("closure_decision", {}).get("real_blind_holdout_executed", False)),
        "statistical_margin_rule_ready": bool(p2648.get("closure_decision", {}).get("can_upgrade_p2647_harness", False)),
        "empirical_holdout_can_replace_source_theorem": False,
    }


def route_matrix(atoms: dict[str, bool]) -> list[dict[str, Any]]:
    routes = [
        {
            "route": "normalization_orbit_beta_one",
            "claim": "set beta=1 by choosing length/UV units",
            "required_atoms": {"normalization_gauge_fixed", "canonical_length_unit_exported"},
            "available_atoms": {"normalization_orbit_exists"},
            "verdict_if_missing": "gauge representative only, not source",
        },
        {
            "route": "information_flux_conservation_beta_one",
            "claim": "derive beta=1 from information-flux conservation",
            "required_atoms": {"target_independent_conservation_constant_exported", "canonical_length_unit_exported", "normalization_gauge_fixed"},
            "available_atoms": {"monotone_flux_functional_exists"},
            "verdict_if_missing": "calibrates to whichever beta supplies the constant",
        },
        {
            "route": "edge_of_chaos_beta_one",
            "claim": "derive beta=1 as unique critical/entropy optimum",
            "required_atoms": {"unique_beta_one_critical_extremum_exported", "canonical_length_unit_exported"},
            "available_atoms": {"neural_attention_prism_available"},
            "verdict_if_missing": "heuristic neural analogy only",
        },
        {
            "route": "micro_zbeta_bridge_source",
            "claim": "derive strict beta=1 or Z_beta=100 from micro renormalization",
            "required_atoms": {"positive_zbeta_source_exported", "bridge_zbeta_source_accepted", "normalization_gauge_fixed"},
            "available_atoms": {"micro_zbeta_measurements_exist"},
            "verdict_if_missing": "micro/strict mismatch and normalization gauge remain",
        },
        {
            "route": "blind_empirical_compression_holdout",
            "claim": "use frozen-kernel empirical compression signature as beta source",
            "required_atoms": {"holdout_real_blind_data_executed", "empirical_holdout_can_replace_source_theorem"},
            "available_atoms": {"holdout_harness_ready", "statistical_margin_rule_ready"},
            "verdict_if_missing": "empirical support route only, never source theorem by itself",
        },
    ]
    rows = []
    for route in routes:
        available_truth = {atom: atoms.get(atom, atom in route["available_atoms"]) for atom in sorted(route["available_atoms"] | route["required_atoms"])}
        missing = sorted(atom for atom in route["required_atoms"] if not available_truth.get(atom, False))
        rows.append({
            "route": route["route"],
            "claim": route["claim"],
            "required_atoms": sorted(route["required_atoms"]),
            "available_atoms": sorted(route["available_atoms"]),
            "atom_truth": available_truth,
            "missing_required_atoms": missing,
            "passes_as_target_independent_beta_source": not missing,
            "verdict": "PASS" if not missing else route["verdict_if_missing"],
        })
    return rows


def closure_decision(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["route"] for row in matrix if row["passes_as_target_independent_beta_source"]]
    return {
        "decision": "NO_TARGET_INDEPENDENT_BETA_SOURCE_ROUTE_PASSES__NORMALIZATION_ORBIT_AND_EMPIRICAL_ROUTES_DO_NOT_EXPORT_SOURCE",
        "passing_beta_source_routes": passing,
        "professorial_verdict": (
            "P2649 separates three things that kept getting conflated: choosing beta=1 as a normalization representative, recovering beta=1 from a chosen tail-ratio target, "
            "and exporting a target-independent beta source theorem.  The first two are algebraically valid but circular as source claims; the empirical P2647/P2648 route is useful support but cannot replace the source theorem."
        ),
        "professorial_closure_path": [
            "First fix a canonical length/UV unit from nadsoliton dynamics, not from the strict kernel target.",
            "Then derive an independent dimensionless conservation constant or operator identity whose unique solution is beta=1 before comparing to K_strict_gate.",
            "Only after that rerun micro Z_beta and role-transfer; do not treat a blind empirical pass as a beta-source substitute.",
            "If no such atom appears, keep beta=1 as a robust working normalization plus empirical-compression hypothesis, not as full ToE sourcehood.",
        ],
        "next_honest_step": (
            "Build a canonical-length/UV-unit theorem or dimensionless conservation identity independent of K_strict_gate; otherwise run the P2647/P2648 holdout only as empirical support while beta-source remains blocked."
        ),
        "beta_source_exported_now": False,
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    decision = payload["closure_decision"]
    lines = [
        "# P2649/S1599 beta-source route decision matrix and normalization-orbit no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps beta-source identity, normalization orbit, micro/Z_beta mismatch, empirical-not-source, and closure guard content before adding the matrix.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Algebraic no-go witnesses",
        "",
        payload["normalization_orbit_witness"]["theorem"],
        payload["normalization_orbit_witness"]["consequence"],
        "",
        payload["tail_ratio_inversion_witness"]["formula"],
        payload["tail_ratio_inversion_witness"]["consequence"],
        "",
        "## Route matrix",
        "",
        "| route | passes as beta source? | missing required atoms | verdict |",
        "| --- | ---: | --- | --- |",
    ])
    for row in payload["beta_source_route_matrix"]:
        lines.append(f"| `{row['route']}` | `{row['passes_as_target_independent_beta_source']}` | `{', '.join(row['missing_required_atoms'])}` | {row['verdict']} |")
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing beta-source routes: `{decision['passing_beta_source_routes']}`.",
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
    atoms = upstream_atoms()
    matrix = route_matrix(atoms)
    decision = closure_decision(matrix)
    payload: dict[str, Any] = {
        "status": "P2649_BETA_SOURCE_ROUTE_MATRIX_NO_TARGET_INDEPENDENT_SOURCE_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {name: sha256_file(path) for name, path in SOURCES.items()},
        "upstream_atoms": atoms,
        "normalization_orbit_witness": normalization_orbit_witness(),
        "tail_ratio_inversion_witness": tail_ratio_inversion_witness(),
        "beta_source_route_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        SOURCES["STRICT_EQUATION_SHEET"],
        "P2649/S1599 beta-source route matrix guard",
        "## P2649/S1599 beta-source route matrix guard\n\n"
        "`P2649/S1599` audits the post-P2648 beta-source routes algebraically.  Under `d'=a*d`, the denominator orbit sends `beta -> beta/a^eta`, so every positive beta can be represented with bare `beta=1` after choosing `a=beta^(1/eta)`; likewise any single tail-ratio target recovers a beta by `beta=(1-q)/(q*b^eta-a^eta)`.  These are normalization/calibration facts, not target-independent source theorems.  The route matrix leaves normalization, flux/criticality, micro `Z_beta`, and empirical holdout routes blocked as beta sources; no bridge completion, role transfer, `QW-2191` discharge, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        SOURCES["STRICT_LAGRANGIAN_DRAFT"],
        "P2649/S1599 beta-source route matrix Ltotal guard",
        "## P2649/S1599 beta-source route matrix Ltotal guard\n\n"
        "`P2649/S1599` does not re-enable `L_total`: beta=1 remains a robust working normalization/compression parameter until a canonical length/UV unit plus target-independent conservation/operator identity is proved; empirical holdout success cannot substitute for that source theorem.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
