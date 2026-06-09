#!/usr/bin/env python3
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
OUT = GEN / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.json"
MD = GEN / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.md"

REPO_ROOT = REPO
QW2064_JSON = REPO_ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json"
QW2064_SCRIPT = REPO_ROOT / "QW_2064_MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE.py"

SOURCE_FILES = {
    "P2625_DAMPING_CLASSIFICATION": GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json",
    "QW2064_MICRO_CONSTANTS_JSON": QW2064_JSON,
    "QW2064_MICRO_CONSTANTS_SCRIPT": QW2064_SCRIPT,
    "STRICT_EQUATION_SHEET": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "LTOTAL_DRAFT": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

NEGATIVE_EXPORT_FLAGS = [
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
    return {"count": len(lines), "samples": lines[:40]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: research terms, not ticket-number lookup.
    patterns = {
        "micro_renormalization_content": "micro-derived|micro derivation|renormalization constants|Z_beta|z_beta|delta_eta|beta_uv",
        "coefficient_source_content": "positive beta|coefficient source|beta renormalization|Z_beta=|beta/beta_tors|strict beta|target.*z_beta",
        "noncircularity_content": "no sector retune|no new optimization|no.*scan|target.*kernel|frozen kernel|selected_kernel|Stage-C",
        "dispersion_interval_content": "CI|IQR|quantile|q25|q50|q75|wide|dispersion|coverage|median",
        "bridge_guard_content": "completion map|bridge-source cut|nonlinear damping completion|role-bearing L_total|QW-2191|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first semantic audit for micro Z_beta source", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def qw2064_metrics(qw2064: dict[str, Any]) -> dict[str, Any]:
    target = float(qw2064["targets"]["z_beta_target"])
    median = float(qw2064["micro_global"]["z_beta_median"])
    q25, q50, q75 = [float(x) for x in qw2064["dispersion_diagnostics"]["z_beta_bin_q25_q50_q75"]]
    log_iqr = float(qw2064["dispersion_diagnostics"]["z_beta_log_iqr"])
    abs_log_ratio = float(qw2064["deviation"]["abs_log_zbeta_ratio"])
    return {
        "target_z_beta": target,
        "micro_global_median_z_beta": median,
        "target_over_median": target / median,
        "median_over_target": median / target,
        "abs_log_median_over_target": abs_log_ratio,
        "z_beta_q25_q50_q75": [q25, q50, q75],
        "target_inside_q25_q75": q25 <= target <= q75,
        "z_beta_log_iqr": log_iqr,
        "wide_ci_warning": bool(qw2064.get("ci_warning", False)),
        "verdict": str(qw2064.get("verdict", "UNKNOWN")),
        "target_definition_depends_on_frozen_kernel": True,
        "exact_target_match": median == target,
    }


def acceptance_table(metrics: dict[str, Any]) -> dict[str, Any]:
    gates = {
        "positive_micro_support": metrics["micro_global_median_z_beta"] > 0.0,
        "target_in_reported_interval": metrics["target_inside_q25_q75"],
        "target_independent_of_selected_kernel": not metrics["target_definition_depends_on_frozen_kernel"],
        "exact_or_bridge_tolerance_theorem": metrics["exact_target_match"],
        "narrow_dispersion_no_wide_ci_warning": not metrics["wide_ci_warning"],
    }
    accepting = all(gates.values())
    return {
        "gates": gates,
        "accepts_positive_beta_renormalization_source": accepting,
        "failed_gates": [name for name, value in gates.items() if not value],
        "interpretation": (
            "QW-2064 gives positive, interval-level support near the requested coefficient, but it is not a strict source theorem: "
            "the target is defined from the selected strict kernel, the reported median is not exactly 100, and the CI is explicitly wide."
        ),
    }


def p2625_update(source_accepts: bool) -> dict[str, Any]:
    return {
        "positive_beta_renormalization_source_after_p2626": source_accepts,
        "nonlinear_damping_completion_source_after_p2626": False,
        "reason": (
            "P2626 audits the best available micro Z_beta lane. Because it does not export the positive coefficient source, "
            "P2625's conditional nonlinear completion schema remains conditional."
        ),
    }


def theorem_export(metrics: dict[str, Any], table: dict[str, Any]) -> dict[str, Any]:
    return {
        "theorem_name": "P2626_T1_micro_zbeta_source_nonpromotion_audit",
        "positive_content": [
            "QW-2064 reports a positive micro-derived Z_beta median and the target 100 lies inside its q25-q75 interval.",
            "QW-2064's gate verdict is pass-with-wide-CI-warning, so it is a useful coefficient-source candidate lane.",
        ],
        "negative_content": [
            "The target Z_beta=100 is computed as beta_selected/beta_uv from the selected frozen strict kernel, so the target is not independent of the kernel being bridged.",
            "The micro median is not exactly 100 and the bin q50 is far from 100; this is support, not an exact theorem.",
            "The reported wide-CI warning blocks promotion to a strict positive_beta_renormalization_source.",
        ],
        "current_verdict": "QW2064_MICRO_ZBETA_IS_SUPPORT_NOT_STRICT_SOURCE",
        "failed_source_gates": table["failed_gates"],
        "recommended_next_honest_step": (
            "Build a noncircular coefficient-source theorem for Z_beta: either derive Z_beta=100 from a target-independent micro operator/normalization law, "
            "or downgrade P2625 to an interval-valued bridge and prove an explicit tolerance theorem before any bridge/role rerun."
        ),
        "metrics_summary": metrics,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    qw2064 = load_json(QW2064_JSON)
    p2625 = load_json(SOURCE_FILES["P2625_DAMPING_CLASSIFICATION"])
    metrics = qw2064_metrics(qw2064)
    table = acceptance_table(metrics)
    theorem = theorem_export(metrics, table)
    payload = {
        "packet_id": "P2626",
        "slice_id": "S1576",
        "status": "P2626_MICRO_ZBETA_SUPPORT_NONPROMOTION_NO_POSITIVE_BETA_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "inherited_artifact_status": {
            "P2625": p2625.get("status", "UNKNOWN"),
            "QW2064": qw2064.get("verdict", "UNKNOWN"),
        },
        "qw2064_metrics": metrics,
        "positive_beta_source_acceptance": table,
        "theorem_export": theorem,
        "p2625_update": p2625_update(table["accepts_positive_beta_renormalization_source"]),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    m = payload["qw2064_metrics"]
    table = payload["positive_beta_source_acceptance"]
    theorem = payload["theorem_export"]
    lines = [
        "# P2626/S1576 Micro Z_beta source nonpromotion audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication grep audit",
        "",
        f"Mode: `{payload['semantic_rg_antiduplication_audit']['mode']}`.",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits; samples retained in JSON certificate.")
    lines.extend([
        "",
        "## QW-2064 coefficient-source audit",
        "",
        f"- target Z_beta: `{m['target_z_beta']}`",
        f"- micro global median Z_beta: `{m['micro_global_median_z_beta']}`",
        f"- median/target: `{m['median_over_target']}`",
        f"- q25/q50/q75: `{m['z_beta_q25_q50_q75']}`",
        f"- target inside q25-q75: `{m['target_inside_q25_q75']}`",
        f"- wide CI warning: `{m['wide_ci_warning']}`",
        f"- target depends on selected kernel: `{m['target_definition_depends_on_frozen_kernel']}`",
        "",
        "## Acceptance verdict",
        "",
        f"Accepts positive_beta_renormalization_source: `{table['accepts_positive_beta_renormalization_source']}`.",
        f"Failed gates: `{table['failed_gates']}`.",
        "",
        "Positive content:",
    ])
    for item in theorem["positive_content"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("Negative content:")
    for item in theorem["negative_content"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Recommended next honest step",
        "",
        theorem["recommended_next_honest_step"],
        "",
        "## Negative export flags",
        "",
    ])
    for flag, value in payload["negative_export_flags"].items():
        lines.append(f"- `{flag}`: `{value}`")
    lines.append("")
    MD.write_text("\n".join(lines), encoding="utf-8")


def update_docs() -> None:
    equation_note = """
## P2626/S1576 micro Z_beta source nonpromotion guard

`P2626/S1576` audits the best current micro-derived coefficient lane for the missing P2625 atom.  QW-2064 gives positive interval-level support for `Z_beta` and includes `100` inside its reported broad interval, but it does not export `positive_beta_renormalization_source`: the target is computed from the already selected strict kernel, the micro median is not exactly `100`, and the report carries a wide-CI warning.  The next admissible step is a noncircular `Z_beta` coefficient-source theorem or an explicit interval/tolerance bridge theorem, not bridge/role promotion.
"""
    ltotal_note = """
## P2626/S1576 micro Z_beta Ltotal guard

`P2626/S1576` keeps role-bearing `L_total` closed.  Micro support for `Z_beta` remains useful evidence, but without a target-independent exact coefficient-source theorem or a separate tolerance bridge theorem it cannot repair P2625, P2620, role transfer, `QW-2191`, or ToE closure.
"""
    append_once(ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md", "P2626/S1576 micro Z_beta source nonpromotion guard", equation_note)
    append_once(ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md", "P2626/S1576 micro Z_beta Ltotal guard", ltotal_note)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    update_docs()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "accepts_positive_beta_source": payload["positive_beta_source_acceptance"]["accepts_positive_beta_renormalization_source"],
        "failed_gates": payload["positive_beta_source_acceptance"]["failed_gates"],
        "recommended_next": payload["theorem_export"]["recommended_next_honest_step"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True, ensure_ascii=False))


if __name__ == "__main__":
    main()
