#!/usr/bin/env python3
"""
QW-2067: Strengthened strict first-principles internal closure gate.

Aggregates:
- QW-2065 strict internal closure pass.
- QW-2066 compatibility-filtered tightening of micro-constants dispersion.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2067_strict_first_principles_internal_closure_strengthened_gate.json"
OUT_MD = ROOT / "RAPORT_QW2067_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_GATE.md"


def main() -> None:
    r2065 = json.loads((ROOT / "report_qw2065_strict_first_principles_internal_closure_gate.json").read_text(encoding="utf-8"))
    r2066 = json.loads((ROOT / "report_qw2066_compatibility_filtered_micro_constants_tightening.json").read_text(encoding="utf-8"))

    flags = {
        "qw2065_strict_internal_pass": bool(r2065.get("strict_internal_pass", False)),
        "qw2066_tightening_pass": bool(r2066.get("verdict") == "COMPATIBILITY_FILTERED_MICRO_CONSTANTS_TIGHTENING_PASS"),
        "qw2066_warning_resolved": bool(r2066.get("tightened_warning_resolved", False)),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    strengthened_pass = bool(all(flags.values()))

    if strengthened_pass:
        verdict = "STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_PASS"
    else:
        verdict = "STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_STRENGTHENED_PARTIAL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "closure_gate": "report_qw2065_strict_first_principles_internal_closure_gate.json",
            "tightening_gate": "report_qw2066_compatibility_filtered_micro_constants_tightening.json",
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "strengthened_pass": strengthened_pass,
        "baseline_ci_warning": bool(r2065.get("foundation_ci_warning", False)),
        "tightened_warning_resolved": bool(r2066.get("tightened_warning_resolved", False)),
        "verdict": verdict,
        "required_next_step": "RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2067: STRICT FIRST-PRINCIPLES INTERNAL CLOSURE STRENGTHENED GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- strengthened_pass: {strengthened_pass}",
        f"- baseline_ci_warning (QW-2065): {out['baseline_ci_warning']}",
        f"- tightened_warning_resolved (QW-2066): {out['tightened_warning_resolved']}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend([
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ])

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2067] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2067] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2067] verdict={verdict} pass_count={pass_count}/{total_flags} strengthened_pass={strengthened_pass}")


if __name__ == "__main__":
    main()
