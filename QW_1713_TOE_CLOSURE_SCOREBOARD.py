#!/usr/bin/env python3
"""
QW-1713: FIN ToE closure scoreboard.

Aggregates new closure-focused studies (QW-1710..1715) and assigns a
transparent readiness status with next experiments.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1713_toe_closure_scoreboard.json"
OUT_MD = ROOT / "RAPORT_QW1713_TOE_CLOSURE_SCOREBOARD.md"


def load_json(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def main() -> None:
    j1710 = load_json(ROOT / "report_qw1710_flavor_topological_operator.json")
    j1711 = load_json(ROOT / "report_qw1711_mass_oos_closure_test.json")
    j1712 = load_json(ROOT / "report_qw1712_ir_lorentz_recovery_test.json")
    j1714 = load_json(ROOT / "report_qw1714_flavor_transition_operator.json")
    j1715 = load_json(ROOT / "report_qw1715_mass_bootstrap_oos.json")
    j1716 = load_json(ROOT / "report_qw1716_uncertainty_budget.json")
    j1717 = load_json(ROOT / "report_qw1717_identifiability_test.json")
    j1718 = load_json(ROOT / "report_qw1718_prospective_oos_benchmark.json")
    j1719 = load_json(ROOT / "report_qw1719_flavor_locality_group_constraint.json")
    j1720 = load_json(ROOT / "report_qw1720_flavor_extended_operator.json")
    j1721 = load_json(ROOT / "report_qw1721_strict_no_leakage_oos.json")
    j1722 = load_json(ROOT / "report_qw1722_parameter_stability_perturbation.json")
    j1723 = load_json(ROOT / "report_qw1723_unified_effective_hamiltonian_integration.json")

    checks = []

    # Flavor closure
    if j1710 is None:
        checks.append(("Flavor CKM/PMNS", "MISSING", 0.0, "Brak raportu QW-1710"))
    else:
        v = j1710.get("verdict", "UNKNOWN")
        if v == "FLAVOR_CLOSURE_DERIVED_PLAUSIBLE":
            checks.append(("Flavor CKM/PMNS", "PASS_DERIVED", 1.0, v))
        elif v == "FLAVOR_CLOSURE_REQUIRES_CALIBRATED_OPERATOR":
            checks.append(("Flavor CKM/PMNS", "PARTIAL_CALIBRATED", 0.5, v))
        else:
            checks.append(("Flavor CKM/PMNS", "FAIL_NOT_CLOSED", 0.0, v))

    # Mass OOS closure
    if j1711 is None:
        checks.append(("Mass OOS", "MISSING", 0.0, "Brak raportu QW-1711"))
    else:
        v = j1711.get("verdict", "UNKNOWN")
        if v == "MASS_SECTOR_OOS_CLOSED":
            checks.append(("Mass OOS", "PASS", 1.0, v))
        else:
            checks.append(("Mass OOS", "FAIL", 0.0, v))

    # Flavor closure with structured operator
    if j1714 is None:
        checks.append(("Flavor Structured Operator", "MISSING", 0.0, "Brak raportu QW-1714"))
    else:
        v = j1714.get("verdict", "UNKNOWN")
        if v == "FLAVOR_CLOSED_WITH_STRUCTURED_OPERATOR":
            checks.append(("Flavor Structured Operator", "PASS", 1.0, v))
        else:
            checks.append(("Flavor Structured Operator", "FAIL", 0.0, v))

    # Mass robustness (bootstrap OOS)
    if j1715 is None:
        checks.append(("Mass Bootstrap OOS", "MISSING", 0.0, "Brak raportu QW-1715"))
    else:
        v = j1715.get("verdict", "UNKNOWN")
        if v == "MASS_OOS_ROBUST":
            checks.append(("Mass Bootstrap OOS", "PASS", 1.0, v))
        else:
            checks.append(("Mass Bootstrap OOS", "FAIL", 0.0, v))

    # IR Lorentz recovery
    if j1712 is None:
        checks.append(("IR Lorentz", "MISSING", 0.0, "Brak raportu QW-1712"))
    else:
        v = j1712.get("verdict", "UNKNOWN")
        if v == "IR_LORENTZ_RECOVERY_OK":
            checks.append(("IR Lorentz", "PASS", 1.0, v))
        else:
            checks.append(("IR Lorentz", "FAIL", 0.0, v))

    # Uncertainty budget readiness
    if j1716 is None:
        checks.append(("Uncertainty Budget", "MISSING", 0.0, "Brak raportu QW-1716"))
    else:
        checks.append(("Uncertainty Budget", "PASS", 1.0, "QW-1716 delivered"))

    # Identifiability
    if j1717 is None:
        checks.append(("Identifiability", "MISSING", 0.0, "Brak raportu QW-1717"))
    else:
        v = j1717.get("verdict", "UNKNOWN")
        if v == "IDENTIFIABLE_WITHIN_TEST_GRID":
            checks.append(("Identifiability", "PASS", 1.0, v))
        elif v == "PARTIALLY_NON_IDENTIFIABLE":
            checks.append(("Identifiability", "PARTIAL", 0.5, v))
        else:
            checks.append(("Identifiability", "FAIL", 0.0, v))

    # Prospective benchmark
    if j1718 is None:
        checks.append(("Prospective OOS Benchmark", "MISSING", 0.0, "Brak raportu QW-1718"))
    else:
        v = j1718.get("verdict", "UNKNOWN")
        if v == "PROSPECTIVE_BENCHMARK_PASS":
            checks.append(("Prospective OOS Benchmark", "PASS", 1.0, v))
        elif v == "PROSPECTIVE_BENCHMARK_PARTIAL":
            checks.append(("Prospective OOS Benchmark", "PARTIAL", 0.5, v))
        else:
            checks.append(("Prospective OOS Benchmark", "FAIL", 0.0, v))

    # Flavor locality + group constraints
    if j1719 is None:
        checks.append(("Flavor Locality+Group", "MISSING", 0.0, "Brak raportu QW-1719"))
    else:
        v = j1719.get("verdict", "UNKNOWN")
        if v == "FLAVOR_LOCALITY_GROUP_CLOSED":
            checks.append(("Flavor Locality+Group", "PASS", 1.0, v))
        else:
            checks.append(("Flavor Locality+Group", "FAIL", 0.0, v))

    # Extended flavor operator
    if j1720 is None:
        checks.append(("Flavor Extended Operator", "MISSING", 0.0, "Brak raportu QW-1720"))
    else:
        v = j1720.get("verdict", "UNKNOWN")
        if v == "FLAVOR_EXTENDED_OPERATOR_CLOSED":
            checks.append(("Flavor Extended Operator", "PASS", 1.0, v))
        else:
            checks.append(("Flavor Extended Operator", "FAIL", 0.0, v))

    # Strict no-leakage OOS
    if j1721 is None:
        checks.append(("Strict OOS No-Leakage", "MISSING", 0.0, "Brak raportu QW-1721"))
    else:
        v = j1721.get("verdict", "UNKNOWN")
        if v == "STRICT_OOS_PASS":
            checks.append(("Strict OOS No-Leakage", "PASS", 1.0, v))
        else:
            checks.append(("Strict OOS No-Leakage", "FAIL", 0.0, v))

    # Parameter stability under perturbations
    if j1722 is None:
        checks.append(("Parameter Stability", "MISSING", 0.0, "Brak raportu QW-1722"))
    else:
        v = j1722.get("verdict", "UNKNOWN")
        if v == "PARAMETERS_STABLE_UNDER_PERTURBATION":
            checks.append(("Parameter Stability", "PASS", 1.0, v))
        else:
            checks.append(("Parameter Stability", "FAIL", 0.0, v))

    # Unified effective Hamiltonian integration
    if j1723 is None:
        checks.append(("Unified Mass+Flavor Hamiltonian", "MISSING", 0.0, "Brak raportu QW-1723"))
    else:
        v = j1723.get("verdict", "UNKNOWN")
        if v == "UNIFIED_EFFECTIVE_HAMILTONIAN_CLOSED":
            checks.append(("Unified Mass+Flavor Hamiltonian", "PASS", 1.0, v))
        else:
            checks.append(("Unified Mass+Flavor Hamiltonian", "FAIL", 0.0, v))

    score = sum(x[2] for x in checks) / max(len(checks), 1)
    if score >= 0.85:
        readiness = "TOE_NEAR_CLOSURE"
    elif score >= 0.50:
        readiness = "TOE_PARTIAL_CLOSURE"
    else:
        readiness = "TOE_OPEN_MAJOR_GAPS"

    next_actions = [
        "QW-1724: Blind external benchmark na nowych danych (bez ponownej kalibracji).",
        "QW-1725: Bayesian model comparison (integrated model vs model rozdzielny).",
        "QW-1726: Symboliczna derivacja mapowania lambda->operator flavor (bez empirycznej ansatzy).",
        "QW-1727: Rozszerzenie testu OOS na wiekszy zestaw observabli i jawny error budget.",
    ]

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": [
            {"domain": c[0], "status": c[1], "score": c[2], "note": c[3]}
            for c in checks
        ],
        "global_score": score,
        "readiness": readiness,
        "next_actions": next_actions,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1713: TOE CLOSURE SCOREBOARD",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {score:.3f}",
        f"- Readiness: **{readiness}**",
        "",
        "## Wyniki domenowe",
    ]
    for c in output["checks"]:
        lines.append(f"- {c['domain']}: {c['status']} ({c['score']:.2f}) | {c['note']}")
    lines.extend(
        [
            "",
            "## Kolejne badania domykające (proponowane)",
        ]
    )
    for a in next_actions:
        lines.append(f"- {a}")
    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1713] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1713] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
