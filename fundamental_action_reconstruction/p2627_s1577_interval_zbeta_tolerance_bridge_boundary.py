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
OUT = GEN / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.json"
MD = GEN / "p2627_s1577_interval_zbeta_tolerance_bridge_boundary.md"

REPO_ROOT = REPO
QW2064_JSON = REPO_ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json"
P2625_JSON = GEN / "p2625_s1575_nonlinear_damping_completion_source_classification.json"
P2626_JSON = GEN / "p2626_s1576_micro_zbeta_source_nonpromotion_audit.json"

BETA_TORS = 0.01
BETA_STRICT = 1.0
ETA_STRICT = 9.0 / 5.0
Z_TARGET = BETA_STRICT / BETA_TORS
DOMAIN_DMAX = 10.0
STRICT_EPSILON = 0.01
PRACTICAL_EPSILON = 0.15

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
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    # Content-first anti-duplication: these are research-concept patterns, not a ticket-number lookup.
    patterns = {
        "interval_bridge_content": (
            "interval-valued bridge|tolerance bridge|bridge tolerance|approximate bridge|"
            "epsilon.*bridge|error envelope|relative denominator error|admissible tolerance"
        ),
        "micro_coefficient_interval_content": (
            "micro-derived|renormalization constants|Z_beta|z_beta|q25|q75|wide CI|"
            "dispersion|median|CI95|confidence interval"
        ),
        "nonlinear_damping_denominator_content": (
            "nonlinear damping|strict denominator|legacy denominator|beta_tors|d\\^eta|"
            "fractal pushforward|positive beta|coefficient source"
        ),
        "nonfit_noncircular_guard_content": (
            "target-independent|selected strict kernel|frozen kernel|no numerical fit|no fit|"
            "no sector retune|noncircular|source theorem"
        ),
        "closure_guard_content": (
            "role-bearing L_total|role transfer|QW-2191|ToE closure|full bridge|"
            "orientation_odd_selector_source|nonlinear_damping_completion_source"
        ),
    }
    return {"tool": "rg", "mode": "content-first semantic audit for interval Z_beta tolerance bridge", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def relative_denominator_error(z_beta: float, d: float) -> float:
    candidate = 1.0 + (z_beta * BETA_TORS) * (d ** ETA_STRICT)
    strict = 1.0 + BETA_STRICT * (d ** ETA_STRICT)
    return abs(candidate - strict) / strict


def sup_error_on_domain(z_beta: float, dmax: float = DOMAIN_DMAX) -> float:
    # For fixed Z, error(d)=|Z/100-1|*d^eta/(1+d^eta), monotone increasing for d>0.
    return relative_denominator_error(z_beta, dmax)


def allowed_z_interval(epsilon: float, dmax: float = DOMAIN_DMAX) -> dict[str, Any]:
    factor = epsilon * (1.0 + dmax ** ETA_STRICT) / (dmax ** ETA_STRICT)
    return {
        "epsilon": epsilon,
        "domain": f"0 < d <= {dmax:g}",
        "lower": Z_TARGET * (1.0 - factor),
        "upper": Z_TARGET * (1.0 + factor),
        "proof": (
            "For the interval candidate D_Z(d)=1+(Z/100)d^(9/5) and strict denominator D_*(d)=1+d^(9/5), "
            "the relative error is |Z/100-1| d^(9/5)/(1+d^(9/5)).  On 0<d<=D this is maximized at D, "
            "so epsilon-admission requires |Z/100-1| <= epsilon*(1+D^(9/5))/D^(9/5)."
        ),
    }


def qw2064_z_values(qw2064: dict[str, Any]) -> dict[str, float]:
    q25, q50, q75 = [float(x) for x in qw2064["dispersion_diagnostics"]["z_beta_bin_q25_q50_q75"]]
    return {
        "target": float(qw2064["targets"]["z_beta_target"]),
        "micro_global_median": float(qw2064["micro_global"]["z_beta_median"]),
        "bin_q25": q25,
        "bin_q50": q50,
        "bin_q75": q75,
    }


def tolerance_certificate(qw2064: dict[str, Any]) -> dict[str, Any]:
    values = qw2064_z_values(qw2064)
    rows = []
    for name, value in values.items():
        rows.append({
            "source": name,
            "z_beta": value,
            "relative_error_at_dmax": sup_error_on_domain(value),
            "passes_strict_1_percent_tolerance": sup_error_on_domain(value) <= STRICT_EPSILON,
            "passes_practical_15_percent_tolerance": sup_error_on_domain(value) <= PRACTICAL_EPSILON,
        })
    interval_width = values["bin_q75"] - values["bin_q25"]
    strict_interval = allowed_z_interval(STRICT_EPSILON)
    practical_interval = allowed_z_interval(PRACTICAL_EPSILON)
    return {
        "theorem_name": "P2627_T1_interval_zbeta_tolerance_boundary",
        "denominators": {
            "candidate": "D_Z(d)=1+(Z_beta/100)*d^(9/5)",
            "strict": "D_*(d)=1+d^(9/5)",
            "legacy_after_fractal_pushforward_without_Z": "D_legacy_push(d)=1+0.01*d^(9/5)",
        },
        "domain": f"0 < d <= {DOMAIN_DMAX:g}",
        "strict_epsilon_interval": strict_interval,
        "practical_epsilon_interval": practical_interval,
        "qw2064_values": values,
        "rows": rows,
        "reported_q25_q75_interval_width": interval_width,
        "reported_interval_subset_of_strict_1_percent_admission": values["bin_q25"] >= strict_interval["lower"] and values["bin_q75"] <= strict_interval["upper"],
        "reported_interval_subset_of_practical_15_percent_admission": values["bin_q25"] >= practical_interval["lower"] and values["bin_q75"] <= practical_interval["upper"],
        "interpretation": (
            "The exact target Z_beta=100 has zero damping-denominator error by definition.  The QW-2064 micro median is close enough "
            "for a loose 15% finite-window denominator tolerance on d<=10, but neither the median nor the reported q25-q75 interval "
            "satisfies a strict 1% bridge tolerance, and the broad q25-q75 interval is far outside any narrow bridge envelope."
        ),
    }


def source_acceptance(cert: dict[str, Any], qw2064: dict[str, Any]) -> dict[str, Any]:
    median_row = next(row for row in cert["rows"] if row["source"] == "micro_global_median")
    gates = {
        "explicit_tolerance_theorem_present": True,
        "target_independent_of_selected_kernel": False,
        "strict_1_percent_median_pass": bool(median_row["passes_strict_1_percent_tolerance"]),
        "strict_1_percent_reported_interval_pass": bool(cert["reported_interval_subset_of_strict_1_percent_admission"]),
        "no_wide_ci_warning": not bool(qw2064.get("ci_warning", False)),
    }
    return {
        "gates": gates,
        "accepts_positive_beta_renormalization_source": all(gates.values()),
        "accepts_interval_valued_support_only": median_row["passes_practical_15_percent_tolerance"] and not all(gates.values()),
        "failed_gates": [name for name, value in gates.items() if not value],
        "status": "INTERVAL_SUPPORT_ONLY_NO_STRICT_SOURCE",
        "reason": (
            "P2627 supplies the missing tolerance formula, but the available micro coefficient remains target-dependent and broad. "
            "At most it supports a future approximate/interval bridge lane; it does not export the exact positive_beta_renormalization_source required by P2625."
        ),
    }


def current_blockers(source_accepts: bool) -> dict[str, Any]:
    blockers = {
        "positive_beta_renormalization_source": source_accepts,
        "nonlinear_damping_completion_source": False,
        "orientation_odd_selector_source": False,
        "role_transfer_theorem": False,
        "qw2191_selector_discharge": False,
        "global_kernel_finality_theorem": False,
    }
    return {
        "blockers": blockers,
        "full_bridge_closed": all(blockers.values()),
        "recommended_next_honest_step": (
            "Do not promote the interval lane to exact bridge closure.  Next, either narrow the micro Z_beta distribution with a "
            "target-independent operator/normalization law until the reported interval lies inside the chosen epsilon envelope, or explicitly "
            "downgrade the bridge program to an approximate effective-kernel theorem with declared epsilon and domain before any role-transfer audit."
        ),
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["tolerance_certificate"]
    acceptance = payload["source_acceptance"]
    lines = [
        "# P2627/S1577 interval Z_beta tolerance bridge boundary",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps research content rather than only ticket names/numbers before adding the tolerance result.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Theorem boundary",
        "",
        cert["strict_epsilon_interval"]["proof"],
        "",
        f"For `epsilon={STRICT_EPSILON}` and `{cert['domain']}`, admitted `Z_beta` must lie in "
        f"`[{cert['strict_epsilon_interval']['lower']:.6f}, {cert['strict_epsilon_interval']['upper']:.6f}]`.",
        f"For a looser diagnostic `epsilon={PRACTICAL_EPSILON}`, admitted `Z_beta` must lie in "
        f"`[{cert['practical_epsilon_interval']['lower']:.6f}, {cert['practical_epsilon_interval']['upper']:.6f}]`.",
        "",
        "## QW-2064 coefficient rows",
        "",
        "| source | Z_beta | sup relative error on d<=10 | pass 1% | pass 15% |",
        "| --- | ---: | ---: | --- | --- |",
    ])
    for row in cert["rows"]:
        lines.append(
            f"| `{row['source']}` | {row['z_beta']:.6f} | {row['relative_error_at_dmax']:.6f} | "
            f"{row['passes_strict_1_percent_tolerance']} | {row['passes_practical_15_percent_tolerance']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        acceptance["reason"],
        "",
        "P2627 therefore does not repair P2625/P2620 and does not reopen role-bearing `L_total`, role transfer, `QW-2191`, or ToE closure.",
        "",
        "## Recommended next honest step",
        "",
        payload["current_blockers_and_recommendation"]["recommended_next_honest_step"],
        "",
        f"Fingerprint: `{payload['fingerprint_sha256']}`",
        "",
    ])
    MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(parents=True, exist_ok=True)
    qw2064 = load_json(QW2064_JSON)
    p2625 = load_json(P2625_JSON)
    p2626 = load_json(P2626_JSON)
    cert = tolerance_certificate(qw2064)
    acceptance = source_acceptance(cert, qw2064)
    blockers = current_blockers(acceptance["accepts_positive_beta_renormalization_source"])
    payload: dict[str, Any] = {
        "status": "P2627_INTERVAL_ZBETA_TOLERANCE_BOUNDARY_NO_STRICT_SOURCE_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "sources": {
            "QW2064_JSON": rel(QW2064_JSON),
            "P2625_JSON": rel(P2625_JSON),
            "P2626_JSON": rel(P2626_JSON),
        },
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "upstream_status": {
            "p2625_status": p2625.get("status"),
            "p2626_status": p2626.get("status"),
            "qw2064_verdict": qw2064.get("verdict"),
            "qw2064_ci_warning": qw2064.get("ci_warning"),
        },
        "constants": {
            "beta_tors": BETA_TORS,
            "beta_strict": BETA_STRICT,
            "eta_strict": ETA_STRICT,
            "z_target": Z_TARGET,
            "domain_dmax": DOMAIN_DMAX,
            "strict_epsilon": STRICT_EPSILON,
            "practical_epsilon": PRACTICAL_EPSILON,
        },
        "tolerance_certificate": cert,
        "source_acceptance": acceptance,
        "current_blockers_and_recommendation": blockers,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["fingerprint_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "fingerprint_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    write_markdown(payload)
    append_once(
        ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
        "## P2627/S1577 interval Z_beta tolerance bridge boundary guard",
        "\n## P2627/S1577 interval Z_beta tolerance bridge boundary guard\n\n"
        "`P2627/S1577` supplies the explicit finite-window tolerance formula for an interval-valued `Z_beta` damping bridge: "
        "for `D_Z(d)=1+(Z_beta/100)d^(9/5)` on `0<d<=10`, the relative denominator error is maximized at `d=10`, so an "
        "`epsilon`-bridge requires `|Z_beta/100-1| <= epsilon*(1+10^(9/5))/10^(9/5)`.  QW-2064's micro median gives only loose interval support, "
        "and its broad q25-q75 interval plus target-dependence and wide-CI warning block promotion to `positive_beta_renormalization_source`.  Thus P2625/P2620 remain unrepaired; no role-transfer, role-bearing `L_total`, `QW-2191`, or ToE closure is reopened.\n",
    )
    append_once(
        ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
        "## P2627/S1577 interval Z_beta Ltotal guard",
        "\n## P2627/S1577 interval Z_beta Ltotal guard\n\n"
        "`P2627/S1577` keeps role-bearing `L_total` closed.  It proves only a tolerance boundary for a possible approximate coefficient bridge, not an exact "
        "target-independent `Z_beta` source theorem.  The micro interval remains too broad for strict bridge repair, so the damping atom, full bridge, role-transfer audit, "
        "`QW-2191` discharge, and ToE closure all remain blocked.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "status": result["status"],
        "source_acceptance": result["source_acceptance"],
        "recommended_next_honest_step": result["current_blockers_and_recommendation"]["recommended_next_honest_step"],
        "out": rel(OUT),
        "md": rel(MD),
    }, indent=2, sort_keys=True, ensure_ascii=False))
