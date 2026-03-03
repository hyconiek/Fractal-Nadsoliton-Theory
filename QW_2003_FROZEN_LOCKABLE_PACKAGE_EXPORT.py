#!/usr/bin/env python3
"""
QW-2003: Export frozen lockable package for true external confirmatory execution.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW2001 = ROOT / "report_qw2001_bounded_gw_triad_lockable_gate.json"
IN_QW2002 = ROOT / "report_qw2002_single_kernel_triple_sector_closure_gate_v3.json"

OUT_PACKAGE = ROOT / "frozen_lockable_triad_package_qw2003.json"
OUT_REPORT_JSON = ROOT / "report_qw2003_frozen_lockable_package_export.json"
OUT_REPORT_MD = ROOT / "RAPORT_QW2003_FROZEN_LOCKABLE_PACKAGE_EXPORT.md"


def canonical_hash(obj: dict) -> str:
    payload = json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    d1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    d2001 = json.loads(IN_QW2001.read_text(encoding="utf-8"))
    d2002 = json.loads(IN_QW2002.read_text(encoding="utf-8"))

    package = {
        "package_id": "QW2003_LOCKABLE_TRIAD_FROZEN_V1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": {
            "q2001_verdict": d2001["verdict"],
            "q2002_verdict": d2002["verdict"],
            "q2002_readiness": d2002["readiness"],
        },
        "frozen_components": {
            "kernel": d2001["frozen_components"]["kernel"],
            "shared_params": d2001["frozen_components"]["shared_params"],
            "gw_operator": d2001["frozen_components"]["gw_operator"],
            "thresholds": d2001["thresholds"],
        },
        "fixed_branch_inputs": {
            "fixed_q_nu": d1969["fixed_q_nu"],
        },
        "frozen_performance_snapshot": {
            "deterministic": d2001["deterministic"],
            "bootstrap": d2001["bootstrap"],
            "local_deterministic_pass_rates": d2001["local_deterministic_pass_rates"],
            "baseline_comparison": d2001["baseline_comparison"],
        },
        "usage_rules": {
            "no_parameter_retuning": True,
            "no_threshold_modification": True,
            "blind_external_only": True,
        },
        "source_reports": [
            IN_QW1969.name,
            IN_QW2001.name,
            IN_QW2002.name,
        ],
    }

    pkg_hash = canonical_hash(package)
    package["package_sha256"] = pkg_hash

    OUT_PACKAGE.write_text(json.dumps(package, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "package_file": OUT_PACKAGE.name,
        "package_sha256": pkg_hash,
        "verdict": "FROZEN_LOCKABLE_PACKAGE_READY",
    }
    OUT_REPORT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2003: FROZEN LOCKABLE PACKAGE EXPORT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Package file: `{OUT_PACKAGE.name}`",
        f"- Package SHA256: `{pkg_hash}`",
        "",
        "## Notes",
        "- package is frozen for external blind confirmatory run,",
        "- parameter retuning is explicitly disallowed.",
        "",
        "## Artifacts",
        f"- JSON report: `{OUT_REPORT_JSON.name}`",
    ]
    OUT_REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2003] Saved package: {OUT_PACKAGE.name}")
    print(f"[QW-2003] Saved JSON:    {OUT_REPORT_JSON.name}")
    print(f"[QW-2003] Saved MD:      {OUT_REPORT_MD.name}")
    print(f"[QW-2003] sha256={pkg_hash}")


if __name__ == "__main__":
    main()
