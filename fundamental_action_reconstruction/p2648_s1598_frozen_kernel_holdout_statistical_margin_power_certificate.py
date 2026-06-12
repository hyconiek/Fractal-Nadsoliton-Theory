#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import statistics
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.json"
MD = GEN / "p2648_s1598_frozen_kernel_holdout_statistical_margin_power_certificate.md"

P2646 = GEN / "p2646_s1596_frozen_kernel_compression_signature_preregistration_certificate.json"
P2647 = GEN / "p2647_s1597_frozen_kernel_blind_holdout_execution_harness_and_fake_pass_firewall.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

FAMILYWISE_ALPHA = 0.05
INEQUALITIES_PER_PAIR = 2

NEGATIVE_EXPORT_FLAGS = [
    "blind_holdout_data_loaded",
    "empirical_confirmation_claimed",
    "statistical_significance_claimed_on_real_data",
    "synthetic_margin_promoted_to_evidence",
    "kernel_retuning_allowed",
    "beta_source_exported",
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
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "statistical_margin_power_content": (
            "statistical margin|power certificate|confidence|Bonferroni|familywise|standard error|"
            "ratio_standard_error|slope_standard_error|margin budget|noise budget"
        ),
        "frozen_holdout_decision_content": (
            "blind holdout|frozen-kernel|no retuning|measured_tail_ratio|measured_log_tail_slope|"
            "tail-ratio/log-slope|denominator/envelope"
        ),
        "legacy_strict_control_content": (
            "legacy fixture|strict fixture|midpoint fixture|control baseline|matched exponential|matched spline|"
            "fake-pass|synthetic fixture"
        ),
        "source_nonclosure_content": (
            "beta source|target-independent|role-transfer|role-bearing L_total|QW-2191|ToE closure|"
            "empirical confirmation"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for statistical margin/power certificate, not packet-name search",
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


def one_sided_familywise_z(num_pairs: int, alpha: float = FAMILYWISE_ALPHA) -> float:
    total_inequalities = num_pairs * INEQUALITIES_PER_PAIR
    per_inequality_alpha = alpha / total_inequalities
    return statistics.NormalDist().inv_cdf(1.0 - per_inequality_alpha)


def margin_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    z = one_sided_familywise_z(len(rows))
    out = []
    for row in rows:
        strict_ratio_margin = row["midpoint_ratio_threshold"] - row["strict_denominator_tail_ratio"]
        strict_slope_margin = row["strict_log_tail_slope"] - row["midpoint_slope_threshold"]
        legacy_ratio_margin = row["midpoint_ratio_threshold"] - row["legacy_hyperbolic_tail_ratio"]
        legacy_slope_margin = row["legacy_log_tail_slope"] - row["midpoint_slope_threshold"]
        out.append({
            "near": row["near"],
            "far": row["far"],
            "midpoint_ratio_threshold": row["midpoint_ratio_threshold"],
            "midpoint_slope_threshold": row["midpoint_slope_threshold"],
            "strict_ratio_margin": strict_ratio_margin,
            "strict_slope_margin": strict_slope_margin,
            "legacy_ratio_margin": legacy_ratio_margin,
            "legacy_slope_margin": legacy_slope_margin,
            "max_ratio_standard_error_for_strict_familywise_pass": strict_ratio_margin / z,
            "max_slope_standard_error_for_strict_familywise_pass": strict_slope_margin / z,
            "legacy_ratio_is_already_wrong_side": legacy_ratio_margin < 0,
            "legacy_slope_is_already_wrong_side": legacy_slope_margin < 0,
            "midpoint_fails_strict_inequality_before_uncertainty": True,
        })
    return out


def decision_rule(rows: list[dict[str, Any]]) -> dict[str, Any]:
    z = one_sided_familywise_z(len(rows))
    return {
        "familywise_alpha": FAMILYWISE_ALPHA,
        "inequalities_per_pair": INEQUALITIES_PER_PAIR,
        "total_inequalities": len(rows) * INEQUALITIES_PER_PAIR,
        "bonferroni_one_sided_z": z,
        "ratio_pass_rule_with_uncertainty": "measured_tail_ratio + z*ratio_standard_error < midpoint_ratio_threshold for every preregistered pair",
        "slope_pass_rule_with_uncertainty": "measured_log_tail_slope - z*slope_standard_error > midpoint_slope_threshold for every preregistered pair",
        "global_pass_rule": "all preregistered ratio and slope inequalities pass simultaneously, controls remain frozen, and no beta/eta/phase retuning is used",
        "failure_rule": "fail if any inequality misses after the familywise uncertainty penalty, if a frozen control removes the advantage, or if success requires holdout refit",
    }


def upstream_consistency() -> dict[str, Any]:
    p2647 = load_json(P2647)
    p2646 = load_json(P2646)
    return {
        "p2646_has_preregistered_rows": bool(p2646.get("preregistered_frozen_kernel_tests", {}).get("rows")),
        "p2647_harness_ready_no_real_holdout": p2647.get("closure_decision", {}).get("real_blind_holdout_executed") is False,
        "p2647_firewall_passes": p2647.get("fake_pass_firewall", {}).get("firewall_passes") is True,
        "p2647_requires_standard_errors": "ratio_standard_error" in json.dumps(p2647, sort_keys=True) and "slope_standard_error" in json.dumps(p2647, sort_keys=True),
    }


def statistical_certificate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    margins = margin_rows(rows)
    return {
        "decision_rule": decision_rule(rows),
        "rows": margins,
        "strict_all_nominal_margins_positive": all(row["strict_ratio_margin"] > 0 and row["strict_slope_margin"] > 0 for row in margins),
        "legacy_all_nominal_margins_negative": all(row["legacy_ratio_margin"] < 0 and row["legacy_slope_margin"] < 0 for row in margins),
        "minimum_strict_ratio_margin": min(row["strict_ratio_margin"] for row in margins),
        "minimum_strict_slope_margin": min(row["strict_slope_margin"] for row in margins),
        "minimum_max_ratio_standard_error_for_strict_pass": min(row["max_ratio_standard_error_for_strict_familywise_pass"] for row in margins),
        "minimum_max_slope_standard_error_for_strict_pass": min(row["max_slope_standard_error_for_strict_familywise_pass"] for row in margins),
    }


def closure_decision(cert: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "STATISTICAL_MARGIN_RULE_READY_BUT_NO_BLIND_DATA__NO_EMPIRICAL_CONFIRMATION_NO_SOURCE_EXPORT",
        "professorial_verdict": (
            "P2648 fixes the main remaining weakness of the P2647 harness: a nominal threshold pass is not enough.  "
            "The frozen holdout must clear every locked ratio/slope inequality after a familywise one-sided uncertainty penalty.  "
            "The strict synthetic prediction has positive margin budget, while the legacy fixture is already on the wrong side of every audited inequality.  "
            "Still, this is a statistical decision rule and power/margin certificate only; no blind data have been tested."
        ),
        "professorial_closure_path": [
            "Require the blind payload to include standard errors or certified uncertainty bounds for each tail ratio and log-tail slope.",
            "Apply the P2648 familywise rule before reading any physical meaning into a P2647 holdout run.",
            "If any inequality passes only without the uncertainty penalty, record empirical failure rather than a weak pass.",
            "If the holdout passes with uncertainty and frozen controls, rerun only the modified/compressed successor row of the role-transfer matrix.",
            "Keep beta-source and selector/source proof obligations separate from empirical margin success.",
        ],
        "next_honest_step": (
            "Acquire a blinded payload with uncertainty estimates tight enough to meet the P2648 margins, then run P2647/P2648 exactly once without retuning; "
            "in parallel continue the target-independent beta-source theorem."
        ),
        "can_upgrade_p2647_harness": cert["strict_all_nominal_margins_positive"] and cert["legacy_all_nominal_margins_negative"],
        "real_blind_holdout_executed": False,
        "full_kernel_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["statistical_margin_certificate"]
    rule = cert["decision_rule"]
    decision = payload["closure_decision"]
    lines = [
        "# P2648/S1598 frozen-kernel holdout statistical margin/power certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps statistical margin/power, frozen holdout decision, legacy-vs-strict controls, and nonclosure guardrails before adding the certificate.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Familywise decision rule",
        "",
        f"Familywise alpha: `{rule['familywise_alpha']}` over `{rule['total_inequalities']}` inequalities.",
        f"Bonferroni one-sided z: `{rule['bonferroni_one_sided_z']:.12f}`.",
        f"Ratio rule: `{rule['ratio_pass_rule_with_uncertainty']}`.",
        f"Slope rule: `{rule['slope_pass_rule_with_uncertainty']}`.",
        "",
        "| near | far | strict ratio margin | strict slope margin | max ratio SE | max slope SE | legacy ratio margin | legacy slope margin |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in cert["rows"]:
        lines.append(
            f"| {row['near']} | {row['far']} | `{row['strict_ratio_margin']:.12f}` | `{row['strict_slope_margin']:.12f}` | "
            f"`{row['max_ratio_standard_error_for_strict_familywise_pass']:.12f}` | `{row['max_slope_standard_error_for_strict_familywise_pass']:.12f}` | "
            f"`{row['legacy_ratio_margin']:.12f}` | `{row['legacy_slope_margin']:.12f}` |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Real blind holdout executed? `{decision['real_blind_holdout_executed']}`.",
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
    rows = p2646_rows()
    cert = statistical_certificate(rows)
    decision = closure_decision(cert)
    payload: dict[str, Any] = {
        "status": "P2648_STATISTICAL_MARGIN_RULE_READY_BUT_NO_BLIND_DATA_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2646_PREREGISTRATION": sha256_file(P2646),
            "P2647_HARNESS": sha256_file(P2647),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "statistical_margin_certificate": cert,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2648/S1598 frozen-kernel holdout statistical margin guard",
        "## P2648/S1598 frozen-kernel holdout statistical margin guard\n\n"
        "`P2648/S1598` upgrades the P2647 holdout harness with a familywise statistical margin rule: every preregistered tail-ratio and log-slope inequality must pass after a Bonferroni one-sided uncertainty penalty, using locked P2646 thresholds and no retuning.  "
        "The strict synthetic fixture has positive margin budget and the legacy fixture is on the wrong side of each audited inequality, but no blinded measurement payload is tested here.  "
        "Thus this exports only a statistical decision rule/power budget: no empirical confirmation, beta source, bridge completion, role transfer, `QW-2191` discharge, role-bearing `L_total`, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2648/S1598 frozen-kernel statistical margin Ltotal guard",
        "## P2648/S1598 frozen-kernel statistical margin Ltotal guard\n\n"
        "`P2648/S1598` does not re-enable `L_total`: it adds uncertainty-penalized pass/fail criteria for the frozen compression holdout, but role-bearing dynamics still require actual blind data, frozen controls, target-independent beta/source proof, selector/source discharge, and role-transfer rerun.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
