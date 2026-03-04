#!/usr/bin/env python3
"""
QW-2071: SM+GR full precision closure gate.

Purpose:
- integrate derivation package status (QW-2069) with radiative program status (QW-2070),
- produce a single strict gate verdict for "full SM+GR precision closure".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2071_sm_gr_full_precision_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW2071_SM_GR_FULL_PRECISION_CLOSURE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2070 = load_json("report_qw2070_full_radiative_program_baseline.json")

    entries = r2069["entries"]
    missing_params: List[str] = sorted([e["id"] for e in entries if e.get("status") == "missing"])
    model_formula_only: List[str] = sorted(
        [e["id"] for e in entries if e.get("strict_level") == "model_formula"]
    )
    strict_derived: List[str] = sorted(
        [e["id"] for e in entries if e.get("strict_level") == "strict_internal_gate" and str(e.get("status", "")).startswith("derived")]
    )

    channels = r2070["channels"]
    missing_radiative_channels: List[str] = sorted([c["id"] for c in channels if c.get("status") == "missing"])
    implemented_radiative_channels: List[str] = sorted([c["id"] for c in channels if str(c.get("status", "")).startswith("implemented")])
    nonclosing_radiative_channels: List[str] = sorted(
        [
            c["id"]
            for c in channels
            if str(c.get("status", "")).startswith("implemented") and not bool(c.get("closure_ready", False))
        ]
    )

    gate_flags = {
        "strict_internal_strengthened_pass": bool(r2069.get("strict_internal_strengthened_pass", False)),
        "full_derivation_package_pass": r2069.get("verdict") == "FULL_SM_GR_DERIVATION_PACKAGE_PASS",
        "radiative_program_pass": r2070.get("verdict") == "FULL_RADIATIVE_PROGRAM_PASS",
        "no_missing_parameters": len(missing_params) == 0,
        "all_radiative_channels_implemented": len(missing_radiative_channels) == 0,
    }
    pass_count = sum(1 for v in gate_flags.values() if v)
    total_flags = len(gate_flags)

    full_precision_pass = all(gate_flags.values())

    if full_precision_pass:
        verdict = "SM_GR_FULL_PRECISION_CLOSURE_PASS"
    elif gate_flags["strict_internal_strengthened_pass"]:
        verdict = "SM_GR_FULL_PRECISION_CLOSURE_PARTIAL_STRONG_INTERNAL"
    else:
        verdict = "SM_GR_FULL_PRECISION_CLOSURE_FAIL"

    required_next_steps = []
    if missing_params:
        required_next_steps.append("Derive all currently missing SM+GR parameters from first principles chain.")
    if missing_radiative_channels:
        required_next_steps.append(
            "Implement currently missing radiative channels: " + ", ".join(missing_radiative_channels) + "."
        )
    if nonclosing_radiative_channels:
        required_next_steps.append(
            "Upgrade non-closing baseline radiative channels to closure-ready precision: "
            + ", ".join(nonclosing_radiative_channels)
            + "."
        )
    if model_formula_only:
        required_next_steps.append("Upgrade model-formula-only entries to strict derivation status.")
    required_next_steps.append("Run truly independent multiteam external confirmation package.")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "derivation_package": "report_qw2069_full_sm_gr_derivation_package.json",
            "radiative_program": "report_qw2070_full_radiative_program_baseline.json",
        },
        "gate_flags": gate_flags,
        "pass_count": int(pass_count),
        "total_flags": int(total_flags),
        "strict_derived_parameters": strict_derived,
        "model_formula_only_parameters": model_formula_only,
        "missing_parameters": missing_params,
        "implemented_radiative_channels": implemented_radiative_channels,
        "missing_radiative_channels": missing_radiative_channels,
        "nonclosing_radiative_channels": nonclosing_radiative_channels,
        "verdict": verdict,
        "required_next_steps": required_next_steps,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2071: SM+GR FULL PRECISION CLOSURE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Gate pass count: {pass_count}/{total_flags}",
        "",
        "## Gate Flags",
    ]
    for k, v in gate_flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Coverage Summary",
            f"- strict-derived parameters: {len(strict_derived)}",
            f"- model-formula-only parameters: {len(model_formula_only)}",
            f"- missing parameters: {len(missing_params)}",
            f"- implemented radiative channels: {len(implemented_radiative_channels)}",
            f"- missing radiative channels: {len(missing_radiative_channels)}",
            f"- implemented but non-closing radiative channels: {len(nonclosing_radiative_channels)}",
            "",
            "## Required Next Steps",
        ]
    )
    for step in required_next_steps:
        lines.append(f"- {step}")

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2071] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2071] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2071] verdict={verdict} pass_count={pass_count}/{total_flags} "
        f"missing_params={len(missing_params)} missing_channels={len(missing_radiative_channels)}"
    )


if __name__ == "__main__":
    main()
