#!/usr/bin/env python3
"""
QW-1718: Prospective out-of-sample benchmark (frozen protocol).

This script defines a frozen benchmark manifest and evaluates precomputed
results from QW-1710..1717 against fixed thresholds.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1718_prospective_oos_benchmark.json"
OUT_MD = ROOT / "RAPORT_QW1718_PROSPECTIVE_OOS_BENCHMARK.md"


def load(path: Path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def main() -> None:
    manifest = {
        "version": "FIN-Prospective-Benchmark-1",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "thresholds": {
            "flavor_shared_mean_rel_pct_max": 15.0,
            "mass_oos_test_mean_pct_max": 10.0,
            "mass_bootstrap_test_median_pct_max": 10.0,
            "mass_bootstrap_gap_mean_pct_max": 7.0,
            "ir_delta_max_mu_le_1e-3": 1e-3,
            "identifiability_mass_cond_max": 1e5,
        },
    }
    manifest_hash = hashlib.sha256(
        json.dumps(manifest, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()

    j1710 = load(ROOT / "report_qw1710_flavor_topological_operator.json")
    j1711 = load(ROOT / "report_qw1711_mass_oos_closure_test.json")
    j1712 = load(ROOT / "report_qw1712_ir_lorentz_recovery_test.json")
    j1715 = load(ROOT / "report_qw1715_mass_bootstrap_oos.json")
    j1717 = load(ROOT / "report_qw1717_identifiability_test.json")

    checks = []

    # Flavor shared test
    if j1710 is not None:
        val = j1710["shared_operator_best"]["ckm_error"]["mean_rel_pct"]
        val2 = j1710["shared_operator_best"]["pmns_error"]["mean_rel_pct"]
        pass_flag = (val <= manifest["thresholds"]["flavor_shared_mean_rel_pct_max"]) and (
            val2 <= manifest["thresholds"]["flavor_shared_mean_rel_pct_max"]
        )
        checks.append(
            {
                "name": "Flavor shared CKM/PMNS",
                "value": {"ckm": val, "pmns": val2},
                "pass": pass_flag,
                "threshold": manifest["thresholds"]["flavor_shared_mean_rel_pct_max"],
            }
        )
    else:
        checks.append({"name": "Flavor shared CKM/PMNS", "value": None, "pass": False, "threshold": None})

    # Mass OOS single split
    if j1711 is not None:
        val = j1711["metrics"]["test_mean_err_pct"]
        checks.append(
            {
                "name": "Mass OOS test_mean",
                "value": val,
                "pass": val <= manifest["thresholds"]["mass_oos_test_mean_pct_max"],
                "threshold": manifest["thresholds"]["mass_oos_test_mean_pct_max"],
            }
        )
    else:
        checks.append({"name": "Mass OOS test_mean", "value": None, "pass": False, "threshold": None})

    # Mass bootstrap
    if j1715 is not None:
        med = j1715["metrics"]["test_median_pct"]
        gap = j1715["metrics"]["gap_mean_pct"]
        pass_flag = (med <= manifest["thresholds"]["mass_bootstrap_test_median_pct_max"]) and (
            gap <= manifest["thresholds"]["mass_bootstrap_gap_mean_pct_max"]
        )
        checks.append(
            {
                "name": "Mass bootstrap robustness",
                "value": {"test_median_pct": med, "gap_mean_pct": gap},
                "pass": pass_flag,
                "threshold": {
                    "test_median_pct": manifest["thresholds"]["mass_bootstrap_test_median_pct_max"],
                    "gap_mean_pct": manifest["thresholds"]["mass_bootstrap_gap_mean_pct_max"],
                },
            }
        )
    else:
        checks.append({"name": "Mass bootstrap robustness", "value": None, "pass": False, "threshold": None})

    # IR Lorentz
    if j1712 is not None:
        val = j1712["baseline"]["delta_ir_max"]
        checks.append(
            {
                "name": "IR Lorentz delta_max",
                "value": val,
                "pass": val <= manifest["thresholds"]["ir_delta_max_mu_le_1e-3"],
                "threshold": manifest["thresholds"]["ir_delta_max_mu_le_1e-3"],
            }
        )
    else:
        checks.append({"name": "IR Lorentz delta_max", "value": None, "pass": False, "threshold": None})

    # Identifiability condition number
    if j1717 is not None:
        val = j1717["mass_identifiability"]["condition_number_xtx"]
        checks.append(
            {
                "name": "Mass identifiability cond(X^T X)",
                "value": val,
                "pass": val <= manifest["thresholds"]["identifiability_mass_cond_max"],
                "threshold": manifest["thresholds"]["identifiability_mass_cond_max"],
            }
        )
    else:
        checks.append({"name": "Mass identifiability cond(X^T X)", "value": None, "pass": False, "threshold": None})

    passed = sum(1 for c in checks if c["pass"])
    total = len(checks)
    score = passed / total if total else 0.0

    if score >= 0.8:
        verdict = "PROSPECTIVE_BENCHMARK_PASS"
    elif score >= 0.5:
        verdict = "PROSPECTIVE_BENCHMARK_PARTIAL"
    else:
        verdict = "PROSPECTIVE_BENCHMARK_FAIL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": manifest,
        "manifest_hash_sha256": manifest_hash,
        "checks": checks,
        "score": score,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1718: PROSPECTIVE OOS BENCHMARK",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Manifest hash: `{manifest_hash}`",
        f"- Score: {score:.3f} ({passed}/{total})",
        f"- Werdykt: **{verdict}**",
        "",
        "## Checki",
    ]
    for c in checks:
        lines.append(f"- {c['name']}: pass={c['pass']}, value={c['value']}, threshold={c['threshold']}")
    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1718] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1718] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
