#!/usr/bin/env python3
"""
QW-1727: TOE empirical closure gate (theory + GW empirical block).

Aggregates outputs from QW-1720..1726 into a single closure verdict.
No legacy files are modified.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1727_toe_empirical_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW1727_TOE_EMPIRICAL_CLOSURE_GATE.md"


def load(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def main() -> None:
    j1720 = load(ROOT / "report_qw1720_flavor_extended_operator.json")
    j1721 = load(ROOT / "report_qw1721_strict_no_leakage_oos.json")
    j1722 = load(ROOT / "report_qw1722_parameter_stability_perturbation.json")
    j1723 = load(ROOT / "report_qw1723_unified_effective_hamiltonian_integration.json")
    j1724 = load(ROOT / "report_qw1724_gw_method_audit.json")
    j1725 = load(ROOT / "report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    j1726 = load(ROOT / "report_qw1726_gw_fin_projection_retest.json")

    checks = []

    def add(domain: str, ok: bool, note: str, partial: bool = False) -> None:
        if ok:
            checks.append({"domain": domain, "status": "PASS", "score": 1.0, "note": note})
        elif partial:
            checks.append({"domain": domain, "status": "PARTIAL", "score": 0.5, "note": note})
        else:
            checks.append({"domain": domain, "status": "FAIL", "score": 0.0, "note": note})

    # Theoretical closure block
    if j1720 is None:
        add("Flavor Extended Operator", False, "Brak QW-1720")
    else:
        add(
            "Flavor Extended Operator",
            j1720.get("verdict") == "FLAVOR_EXTENDED_OPERATOR_CLOSED",
            j1720.get("verdict", "UNKNOWN"),
        )

    if j1721 is None:
        add("Mass Strict OOS", False, "Brak QW-1721")
    else:
        add("Mass Strict OOS", j1721.get("verdict") == "STRICT_OOS_PASS", j1721.get("verdict", "UNKNOWN"))

    if j1722 is None:
        add("Parameter Stability", False, "Brak QW-1722")
    else:
        add(
            "Parameter Stability",
            j1722.get("verdict") == "PARAMETERS_STABLE_UNDER_PERTURBATION",
            j1722.get("verdict", "UNKNOWN"),
        )

    if j1723 is None:
        add("Unified Effective Hamiltonian", False, "Brak QW-1723")
    else:
        add(
            "Unified Effective Hamiltonian",
            j1723.get("verdict") == "UNIFIED_EFFECTIVE_HAMILTONIAN_CLOSED",
            j1723.get("verdict", "UNKNOWN"),
        )

    # Empirical GW block
    if j1724 is None:
        add("GW Method Audit", False, "Brak QW-1724")
    else:
        risk = j1724.get("total_risk_points")
        ok = j1724.get("verdict") == "GW_PIPELINE_METHOD_HIGH_CONFIDENCE"
        partial = j1724.get("verdict") == "GW_PIPELINE_METHOD_MODERATE_RISK"
        add("GW Method Audit", ok, f"risk_points={risk}", partial=partial)

    if j1725 is None:
        add("GW Strict Reanalysis", False, "Brak QW-1725")
    else:
        add(
            "GW Strict Reanalysis",
            j1725.get("verdict") == "GW_CROSS_HURST_ANOMALY_ROBUST",
            j1725.get("verdict", "UNKNOWN"),
        )

    if j1726 is None:
        add("GW FIN Projection", False, "Brak QW-1726")
    else:
        add(
            "GW FIN Projection",
            j1726.get("verdict") == "FIN_023_TO_031_PROJECTION_SUPPORTED",
            j1726.get("verdict", "UNKNOWN"),
        )

    score = sum(c["score"] for c in checks) / max(len(checks), 1)
    if score >= 0.85:
        readiness = "TOE_CLOSED_WITH_EMPIRICAL_SUPPORT"
    elif score >= 0.60:
        readiness = "TOE_PARTIAL_NEEDS_TARGETED_REPAIRS"
    else:
        readiness = "TOE_OPEN_NOT_EMPIRICALLY_CLOSED"

    # Hard gate: require all core domains
    core_domains = {
        "Flavor Extended Operator": False,
        "Mass Strict OOS": False,
        "Parameter Stability": False,
        "Unified Effective Hamiltonian": False,
        "GW Strict Reanalysis": False,
        "GW FIN Projection": False,
    }
    for c in checks:
        if c["domain"] in core_domains and c["status"] == "PASS":
            core_domains[c["domain"]] = True
    hard_gate_pass = all(core_domains.values())

    if hard_gate_pass:
        final_verdict = "TOE_HARD_GATE_PASS"
    else:
        final_verdict = "TOE_HARD_GATE_FAIL"

    next_actions = [
        "QW-1728: Latency-aware cross-detector analysis (explicit 10ms propagation model in estimator).",
        "QW-1729: Blind holdout epoch benchmark on untouched GPS windows.",
        "QW-1730: Unified symbolic derivation replacing fitted mapping lambda->flavor operator.",
        "QW-1731: External replication package (single-command reproducibility, frozen seeds, frozen manifests).",
    ]

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": score,
        "readiness": readiness,
        "hard_gate_core_pass_flags": core_domains,
        "hard_gate_verdict": final_verdict,
        "next_actions": next_actions,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1727: TOE EMPIRICAL CLOSURE GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {score:.3f}",
        f"- Readiness: **{readiness}**",
        f"- Hard gate: **{final_verdict}**",
        "",
        "## Wyniki domenowe",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} ({c['score']:.2f}) | {c['note']}")
    lines.extend(["", "## Core gate flags"])
    for k, v in core_domains.items():
        lines.append(f"- {k}: {v}")
    lines.extend(["", "## Kolejne badania", *[f"- {x}" for x in next_actions]])
    lines.extend(["", "## Artefakty", f"- JSON szczegolowy: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1727] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1727] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
