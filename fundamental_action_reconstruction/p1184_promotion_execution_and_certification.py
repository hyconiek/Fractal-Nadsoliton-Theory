#!/usr/bin/env python3
"""P1184 promotion execution: run certified strict-E2E on promoted candidate."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1183 = json.loads((GEN / "p1183_candidate_promotion_gate_summary.json").read_text(encoding="utf-8"))
    promoted = p1183.get("winner_candidate")
    allowed = bool(p1183.get("promote_candidate")) and isinstance(promoted, str)

    reg_path = GEN / "p1184_promoted_candidate_registry.json"
    if allowed:
        reg_path.write_text(json.dumps({"candidates": [promoted]}, indent=2) + "\n", encoding="utf-8")

    run = None
    integration = None
    if allowed:
        cmd = [sys.executable, str(ROOT / "p1169_e2e_candidate_integration.py"), "--strict-e2e", "--require-out-of-locality-robustness", "--robustness-threshold", "0.6", "--candidate", promoted, "--registry", str(reg_path)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        run = {"cmd": cmd, "returncode": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}
        integration = json.loads((GEN / "p1169_e2e_candidate_integration_summary.json").read_text(encoding="utf-8"))

    certified = bool(allowed and integration and integration.get("integrated_pass"))
    out = {
        "packet": "P1184",
        "as_of": "2026-05-10",
        "promotion_gate_pass": allowed,
        "promoted_candidate": promoted,
        "promoted_registry": str(reg_path) if allowed else None,
        "integration_run": run,
        "integration_summary": {
            "integrated_pass": integration.get("integrated_pass") if integration else None,
            "robust_fraction": ((integration or {}).get("out_of_locality_robustness_summary") or {}).get("robust_fraction") if integration else None,
            "robustness_threshold": (integration or {}).get("robustness_threshold") if integration else None,
        },
        "certified_for_next_stage": certified,
        "note": "Promotion execution only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1184_promotion_execution_and_certification_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1184] certified_for_next_stage={certified} wrote {out_path}")

if __name__ == "__main__":
    main()
