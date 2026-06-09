#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.json"
MD = GEN / "p2629_s1579_zbeta_normalization_gauge_obstruction.md"

REPO_ROOT = REPO
QW2048_JSON = REPO_ROOT / "report_qw2048_spectral_phase_locked_pointwise_derivation.json"
QW2064_JSON = REPO_ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json"
P2627_JSON = GEN / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.json"
P2628_JSON = GEN / "p2628_s1578_target_blind_micro_zbeta_filter_narrowing_audit.json"

STRICT_BETA = 1.0
UV_BETA = 0.01
STRICT_Z_TARGET = STRICT_BETA / UV_BETA
EXACT_TOL = 0.0
STRICT_REL_TOL = 0.01

NEGATIVE_EXPORT_FLAGS = [
    "uv_normalization_source_exported",
    "positive_beta_renormalization_source_exported",
    "nonlinear_damping_completion_source_exported",
    "full_legacy_to_strict_bridge_revalidated",
    "orientation_odd_selector_source_exported",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.lean", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO_ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: research content, not packet IDs.
    patterns = {
        "normalization_gauge_content": (
            "normalization law|UV normalization|beta_uv|normalization gauge|scale convention|"
            "unit convention|canonical normalization|target-independent normalization"
        ),
        "micro_operator_source_content": (
            "micro operator|microscopic operator|conservation.*identity|normalization identity|"
            "source theorem|coefficient source|positive beta|renormalization source"
        ),
        "zbeta_exactness_content": (
            "Z_beta|z_beta|beta_median|beta_strict|strict beta|beta/beta_tors|"
            "micro-derived renormalization|wide CI|median"
        ),
        "nonfit_target_independence_content": (
            "target-independent|target blind|target-blind|post-hoc|no numerical fit|no fit|"
            "selected strict kernel|frozen kernel|retune"
        ),
        "bridge_closure_guard_content": (
            "nonlinear_damping_completion_source|role-bearing L_total|role transfer|QW-2191|"
            "ToE closure|full bridge|K_legacy_ont|K_strict_gate"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for Z_beta normalization-gauge obstruction", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def qw_metrics(qw2048: dict[str, Any], qw2064: dict[str, Any]) -> dict[str, float]:
    beta_median = float(qw2048["global_estimates"]["beta_median"])
    z_median = float(qw2064["micro_global"]["z_beta_median"])
    z_target = float(qw2064["targets"]["z_beta_target"])
    return {
        "beta_uv": float(qw2064["uv_canonical_constants"]["beta_uv"]),
        "strict_beta": float(qw2064["frozen_kernel"]["beta"]),
        "micro_beta_median": beta_median,
        "z_beta_median": z_median,
        "z_beta_target": z_target,
        "micro_over_strict_beta": beta_median / float(qw2064["frozen_kernel"]["beta"]),
        "z_median_over_z_target": z_median / z_target,
        "abs_beta_error": abs(beta_median - float(qw2064["frozen_kernel"]["beta"])),
        "relative_beta_error": abs(beta_median / float(qw2064["frozen_kernel"]["beta"]) - 1.0),
    }


def normalization_orbit(metrics: dict[str, float]) -> dict[str, Any]:
    # lambda rescales the UV unit beta_uv -> lambda*beta_uv.  Absolute Z values move, but the mismatch ratio is invariant.
    lambdas = [0.5, 1.0, metrics["z_beta_median"] / STRICT_Z_TARGET, 2.0, 10.0]
    rows = []
    for lam in lambdas:
        transformed_uv = lam * metrics["beta_uv"]
        transformed_z_median = metrics["micro_beta_median"] / transformed_uv
        transformed_z_target = metrics["strict_beta"] / transformed_uv
        rows.append({
            "lambda": lam,
            "beta_uv_prime": transformed_uv,
            "z_beta_median_prime": transformed_z_median,
            "z_beta_target_prime": transformed_z_target,
            "median_over_target_invariant": transformed_z_median / transformed_z_target,
            "absolute_target_is_100": transformed_z_target == STRICT_Z_TARGET,
        })
    return {
        "group_action": "R_+ normalization-gauge action: beta_uv -> lambda*beta_uv, Z -> Z/lambda",
        "invariant": "Z_median/Z_target = beta_median/beta_strict",
        "rows": rows,
        "invariant_ratio": metrics["z_median_over_z_target"],
        "proof": (
            "A target-independent source cannot be the bare number Z_beta=100 unless an independent UV/legacy normalization theorem fixes beta_uv=0.01. "
            "Under beta_uv -> lambda beta_uv, both Z_median and Z_target scale by 1/lambda, while their ratio beta_median/beta_strict is invariant. "
            "The current invariant ratio is not 1, so no normalization-gauge choice both proves the absolute value 100 and removes the micro/strict mismatch."
        ),
    }


def exact_source_gate(metrics: dict[str, float], orbit: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "uv_normalization_fixed_by_independent_theorem": False,
        "micro_beta_equals_strict_beta_exactly": metrics["abs_beta_error"] <= EXACT_TOL,
        "micro_beta_within_strict_1pct": metrics["relative_beta_error"] <= STRICT_REL_TOL,
        "normalization_invariant_mismatch_removed": abs(orbit["invariant_ratio"] - 1.0) <= EXACT_TOL,
        "wide_ci_warning_absent": False,
    }
    return {
        "gates": gates,
        "accepts_uv_normalization_source": gates["uv_normalization_fixed_by_independent_theorem"],
        "accepts_positive_beta_renormalization_source": all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "NO_UV_NORMALIZATION_SOURCE_NO_EXACT_POSITIVE_BETA_SOURCE",
        "interpretation": (
            "The remaining damping coefficient problem is not solved by another filter over existing micro rows.  The absolute number Z_beta=100 is tied to the "
            "UV/legacy normalization beta_uv=beta_tors=0.01, while the normalization-invariant content of QW-2048/QW-2064 is the mismatch ratio "
            "beta_median/beta_strict≈1.1474.  A real source theorem must derive beta=1 or beta_uv=0.01 from a target-independent micro identity."
        ),
    }


def recommendation() -> str:
    return (
        "Next honest step: formulate and test a target-independent micro normalization identity for beta itself (for example a conserved flux/action "
        "normalization whose stationary value is beta=1) before mentioning Z_beta=100.  If no such identity is available, stop the exact-bridge lane and "
        "downgrade the damping side to an explicitly approximate effective bridge with declared epsilon/domain; do not rerun role-transfer or L_total."
    )


def write_markdown(payload: dict[str, Any]) -> None:
    metrics = payload["qw_metrics"]
    gate = payload["exact_source_gate"]
    lines = [
        "# P2629/S1579 Z_beta normalization-gauge obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps research content around normalization, source theorems, target independence, and bridge guards rather than only names/numbers.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Normalization-gauge theorem",
        "",
        payload["normalization_orbit_certificate"]["proof"],
        "",
        "## Current invariant numbers",
        "",
        f"- `micro_beta_median = {metrics['micro_beta_median']:.12f}`",
        f"- `strict_beta = {metrics['strict_beta']:.12f}`",
        f"- `Z_beta_median / Z_beta_target = beta_median / beta_strict = {metrics['z_median_over_z_target']:.12f}`",
        f"- relative mismatch to exact strict beta: `{metrics['relative_beta_error']:.12f}`",
        "",
        "## Verdict",
        "",
        gate["interpretation"],
        "",
        f"Failed gates: `{gate['failed_gates']}`.",
        "",
        "P2629 therefore does not export `positive_beta_renormalization_source` and does not repair P2625/P2620, role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure.",
        "",
        "## Recommended next honest step",
        "",
        payload["recommended_next_honest_step"],
        "",
        f"Fingerprint: `{payload['fingerprint_sha256']}`",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(parents=True, exist_ok=True)
    qw2048 = load_json(QW2048_JSON)
    qw2064 = load_json(QW2064_JSON)
    p2627 = load_json(P2627_JSON)
    p2628 = load_json(P2628_JSON)
    metrics = qw_metrics(qw2048, qw2064)
    orbit = normalization_orbit(metrics)
    gate = exact_source_gate(metrics, orbit)
    payload: dict[str, Any] = {
        "status": "P2629_ZBETA_NORMALIZATION_GAUGE_OBSTRUCTION_NO_STRICT_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "sources": {
            "QW2048_JSON": rel(QW2048_JSON),
            "QW2064_JSON": rel(QW2064_JSON),
            "P2627_JSON": rel(P2627_JSON),
            "P2628_JSON": rel(P2628_JSON),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_status": {
            "qw2048_verdict": qw2048.get("verdict"),
            "qw2064_verdict": qw2064.get("verdict"),
            "qw2064_ci_warning": qw2064.get("ci_warning"),
            "p2627_status": p2627.get("status"),
            "p2628_status": p2628.get("status"),
        },
        "qw_metrics": metrics,
        "normalization_orbit_certificate": orbit,
        "exact_source_gate": gate,
        "recommended_next_honest_step": recommendation(),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2629/S1579 Z_beta normalization-gauge obstruction guard",
        "\n## P2629/S1579 Z_beta normalization-gauge obstruction guard\n\n"
        "`P2629/S1579` separates absolute `Z_beta=100` from normalization-invariant micro content.  Under `beta_uv -> lambda beta_uv`, both `Z_beta_median` and `Z_beta_target` scale by `1/lambda`; only `Z_beta_median/Z_beta_target = beta_median/beta_strict` is invariant, and the current invariant ratio is about `1.1474`, not `1`.  Therefore neither UV normalization convention nor target-blind filtering exports `positive_beta_renormalization_source`; P2625/P2620, role-transfer, role-bearing `L_total`, `QW-2191`, and ToE closure remain blocked pending a genuine target-independent micro normalization identity for `beta` itself.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2629/S1579 Z_beta normalization-gauge Ltotal guard",
        "\n## P2629/S1579 Z_beta normalization-gauge Ltotal guard\n\n"
        "`P2629/S1579` keeps role-bearing `L_total` closed.  The exact damping coefficient cannot be obtained by treating the absolute number `Z_beta=100` as source data; the normalization-invariant micro mismatch remains nonzero, so a target-independent beta-normalization theorem is still required before damping completion, role transfer, `QW-2191`, or ToE closure can be rerun.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "invariant_ratio": result["normalization_orbit_certificate"]["invariant_ratio"],
        "failed_gates": result["exact_source_gate"]["failed_gates"],
        "recommended_next": result["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
