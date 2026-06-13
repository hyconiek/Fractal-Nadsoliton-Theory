#!/usr/bin/env python3
"""P2692/S1642: target-independent positive beta/Z_beta source audit.

This executes the post-P2691 recommendation: audit whether the remaining P2680
non-selector damping/compression atom supplies a target-independent positive
beta or Z_beta source, without canonical UV-unit replay, beta_tors->chi11,
selector replay, role transfer, generic bridge completion, L_total, or ToE
promotion.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.json"
MD = GEN / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2691": GEN / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
    "P2629": GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json",
    "P2630": GEN / "p2630_s1580_strict_beta_source_vs_bridge_zbeta_separation.json",
    "P2649": GEN / "p2649_s1599_beta_source_route_decision_matrix_and_normalization_orbit_no_go.json",
    "P2651": GEN / "p2651_s1601_beta_one_gauge_fixed_working_normalization_contract.json",
    "P2653": GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "target_independent_positive_beta_source_exported",
    "target_independent_zbeta_source_exported",
    "canonical_uv_unit_replayed_as_source",
    "beta_tors_chi11_imported",
    "selector_replay_imported",
    "role_transfer_started",
    "generic_bridge_completion_claimed",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "p2691_selected_p2692": r"P2692|target-independent positive beta|positive beta/Z_beta|Z_beta source",
        "beta_zbeta_obstructions": r"P2629|P2630|Z_beta|z_beta|normalization gauge|source separation",
        "normalization_orbit_contract": r"P2649|P2651|beta=1|normalization orbit|gauge-fixed working normalization",
        "typed_metric_uv_obligations": r"P2653|typed nadsoliton metric|UV source|scale-orbit quotient|dimensionless conservation",
        "forbidden_imports": r"canonical UV-unit replay|beta_tors.*chi_?11|selector replay|role transfer|generic bridge|bridge completion|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first target-independent positive beta/Z_beta source audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    p2691 = load_json(INPUTS["P2691"])
    p2680 = load_json(INPUTS["P2680"])
    p2649 = load_json(INPUTS["P2649"])
    p2651 = load_json(INPUTS["P2651"])
    p2653 = load_json(INPUTS["P2653"])
    atoms = {row.get("atom"): row for row in p2680.get("source_atom_inventory", [])}
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2691_selected_p2692": "P2692" in p2691.get("decision", {}).get("next_honest_step", ""),
        "p2680_positive_beta_zbeta_source_exported": atoms.get("target_independent_positive_beta_or_z_beta_source", {}).get("source_theorem_exported") is True,
        "p2680_canonical_length_uv_source_exported": atoms.get("canonical_length_or_uv_unit_source", {}).get("source_theorem_exported") is True,
        "p2649_beta_source_exported": p2649.get("closure_decision", {}).get("beta_source_exported_now") is True,
        "p2649_passing_routes": p2649.get("closure_decision", {}).get("passing_beta_source_routes", []),
        "p2651_beta_one_working_gauge_allowed": p2651.get("closure_decision", {}).get("beta_one_working_gauge_allowed") is True,
        "p2651_beta_source_exported": p2651.get("closure_decision", {}).get("beta_source_exported_now") is True,
        "p2653_scale_orbit_equivalence_verified": p2653.get("closure_decision", {}).get("scale_orbit_equivalence_verified") is True,
        "p2653_beta_source_exported": p2653.get("closure_decision", {}).get("beta_source_exported_now") is True,
        "p2653_canonical_unit_exported": p2653.get("closure_decision", {}).get("canonical_unit_exported_now") is True,
    }


def beta_orbit_witness() -> dict[str, Any]:
    eta = 9.0 / 5.0
    betas = [0.01, 0.0927, 0.92, 1.0, 2.5, 100.0]
    rows = []
    for beta in betas:
        a_to_beta_one = beta ** (1.0 / eta)
        beta_prime = beta / (a_to_beta_one ** eta)
        invariant_samples = []
        for d in [1.0, 2.0, 7.0, 12.0]:
            transformed_d = a_to_beta_one * d
            original = beta * (d ** eta)
            transformed = beta_prime * (transformed_d ** eta)
            invariant_samples.append({"d": d, "original_beta_d_eta": original, "transformed_beta_d_eta": transformed, "residual": transformed - original})
        rows.append({"beta": beta, "a_to_beta_one": a_to_beta_one, "beta_prime": beta_prime, "max_abs_invariant_residual": max(abs(x["residual"]) for x in invariant_samples), "samples": invariant_samples})
    return {
        "eta": eta,
        "formula": "d' = a*d, beta' = beta/a^eta; choosing a=beta^(1/eta) gives beta'=1 for every beta>0",
        "rows": rows,
        "all_positive_betas_have_beta_one_representative": all(abs(row["beta_prime"] - 1.0) < 1e-12 for row in rows),
        "max_abs_invariant_residual": max(row["max_abs_invariant_residual"] for row in rows),
        "source_generated_by_orbit": False,
    }


def tail_ratio_inversion_witness() -> dict[str, Any]:
    eta = 9.0 / 5.0
    a, b = 1.0, 7.0
    qs = [0.05, 0.1, 0.2, 0.5]
    rows = []
    for q in qs:
        beta = (1.0 - q) / (q * (b ** eta) - (a ** eta))
        rows.append({"target_tail_ratio_q": q, "d_pair": [a, b], "recovered_beta": beta, "positive_beta": beta > 0})
    return {
        "formula": "For ratio q=(1+beta*a^eta)/(1+beta*b^eta), beta=(1-q)/(q*b^eta-a^eta)",
        "rows": rows,
        "positive_beta_recoverable_for_multiple_declared_targets": all(row["positive_beta"] for row in rows),
        "target_independent_source_generated": False,
    }


def source_candidate_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "normalization_orbit_beta_equals_one_representative",
            "positive_beta_available": True,
            "target_independent": False,
            "source_exported_now": False,
            "reason": "The orbit proves representability of every positive beta by beta=1 after unit choice; it does not choose the unit/source.",
        },
        {
            "candidate": "micro_Z_beta_renormalization_interpretation",
            "positive_beta_available": True,
            "target_independent": False,
            "source_exported_now": False,
            "reason": "P2629/P2630 separate Z_beta normalization/bridge bookkeeping from a strict beta source theorem.",
        },
        {
            "candidate": "tail_ratio_or_empirical_target_inversion",
            "positive_beta_available": True,
            "target_independent": False,
            "source_exported_now": False,
            "reason": "P2649-style inversion recovers beta only after a declared target or holdout unit map.",
        },
        {
            "candidate": "canonical_length_or_uv_unit_reuse",
            "positive_beta_available": reads["p2680_canonical_length_uv_source_exported"],
            "target_independent": reads["p2653_canonical_unit_exported"],
            "source_exported_now": False,
            "reason": "P2689/P2690 and P2653 leave canonical unit/UV source unexported; replay is forbidden here.",
        },
        {
            "candidate": "dimensionless_conservation_or_operator_identity_unique_beta_one",
            "positive_beta_available": False,
            "target_independent": True,
            "source_exported_now": False,
            "reason": "P2653 names this as an obligation but no current artifact exports the identity.",
        },
    ]


def decision(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row for row in matrix if row["positive_beta_available"] and row["target_independent"] and row["source_exported_now"]]
    return {
        "decision": "P2692_TARGET_INDEPENDENT_POSITIVE_BETA_ZBETA_SOURCE_AUDIT_NO_EXPORTED_SOURCE_NO_FALSE_PASS",
        "passing_source_candidates": [row["candidate"] for row in passing],
        "target_independent_positive_beta_or_zbeta_source_exported_now": bool(passing),
        "bounded_no_go_for_current_beta_zbeta_atom": not passing,
        "professorial_verdict": (
            "P2692 gives the beta/Z_beta atom its strongest finite audit.  The normalization orbit is real and positive: every beta>0 has a beta=1 representative after a distance-unit rescaling, and tail-ratio equations recover positive beta after a declared target.  But both facts are source-insufficient.  Current artifacts keep Z_beta as normalization/bridge bookkeeping, beta=1 as a gauge-fixed working representative, empirical inversion as target-dependent, and canonical UV-unit/conservation identity as unexported.  Thus no target-independent positive beta/Z_beta source is exported on current evidence."
        ),
        "next_honest_step": (
            "P2693 should be a post-P2680 non-selector source-inventory closure/state-map reconciliation: mark amplitude, canonical UV-unit, and positive beta/Z_beta atoms as bounded no-go on current artifacts, then choose any further move only from a fresh state-map object rather than replaying generic bridge completion."
        ),
        "role_transfer_started_now": False,
        "ltotal_promoted_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2692/S1642 target-independent positive beta/Z_beta source audit", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    orbit = payload["beta_orbit_witness"]
    lines.extend(["", "## Beta orbit witness", f"`{orbit['formula']}`", f"all_positive_betas_have_beta_one_representative = `{orbit['all_positive_betas_have_beta_one_representative']}`; max_abs_invariant_residual = `{orbit['max_abs_invariant_residual']}`."])
    inv = payload["tail_ratio_inversion_witness"]
    lines.extend(["", "## Tail-ratio inversion witness", f"`{inv['formula']}`", f"positive_beta_recoverable_for_multiple_declared_targets = `{inv['positive_beta_recoverable_for_multiple_declared_targets']}`; target_independent_source_generated = `{inv['target_independent_source_generated']}`."])
    lines.extend(["", "## Source candidate matrix"])
    for row in payload["source_candidate_matrix"]:
        lines.append(f"- `{row['candidate']}`: positive=`{row['positive_beta_available']}`, target_independent=`{row['target_independent']}`, exported=`{row['source_exported_now']}` — {row['reason']}")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    orbit = beta_orbit_witness()
    inversion = tail_ratio_inversion_witness()
    matrix = source_candidate_matrix(reads)
    payload: dict[str, Any] = {
        "status": "P2692_TARGET_INDEPENDENT_POSITIVE_BETA_ZBETA_SOURCE_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "beta_orbit_witness": orbit,
        "tail_ratio_inversion_witness": inversion,
        "source_candidate_matrix": matrix,
        "decision": decision(matrix),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2692/S1642 target-independent positive beta/Z_beta source audit",
        "## P2692/S1642 target-independent positive beta/Z_beta source audit\n\n"
        "`P2692/S1642` audits the remaining P2680 damping/compression source atom.  The finite orbit calculation confirms that every positive `beta` has a `beta=1` representative under `d' = a*d`, `beta' = beta/a^eta`, and tail-ratio equations can recover positive beta after a declared target.  These are normalization/target facts, not target-independent source theorems: `Z_beta` remains bridge bookkeeping, `beta=1` remains a gauge-fixed working representative, empirical inversion remains target-dependent, and no canonical UV unit or dimensionless conservation/operator identity is exported.  Therefore no positive `beta/Z_beta` source, bridge completion, role transfer, selector closure, role-bearing `L_total`, or ToE closure is claimed.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2692/S1642 beta/Z_beta Ltotal guard",
        "## P2692/S1642 beta/Z_beta Ltotal guard\n\n"
        "`P2692/S1642` keeps `L_total` nonpromoted.  Positivity and gauge representability of `beta` do not create a variational damping coefficient source; the missing typed UV unit or independent operator/conservation identity remains unexported.\n",
    )
    append_once(
        AGENTS,
        "Current beta/Z_beta source guardrail (P2692/S1642, 2026-06-13)",
        "## Current beta/Z_beta source guardrail (P2692/S1642, 2026-06-13)\n\n"
        "- P2692 confirms positive-beta orbit/gauge representability and target-dependent tail-ratio inversion, but finds no target-independent positive `beta/Z_beta` source theorem.\n"
        "- Freeze the current P2680 damping/compression beta-source atom as bounded no-go; the next move should be a post-P2680 non-selector source-inventory closure/state-map reconciliation, not generic bridge replay, role transfer, selector replay, canonical UV-unit replay, `L_total`, or ToE closure.\n",
    )
    return payload


if __name__ == "__main__":
    main()
