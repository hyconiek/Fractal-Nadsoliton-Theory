#!/usr/bin/env python3
"""
QW-2065: Strict first-principles internal closure gate.

Aggregates:
- QW-2063 physical triad closure with deterministic no-scan map.
- QW-2064 foundational renormalization constants support from micro derivation.

No new fitting/scanning is introduced.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2065_strict_first_principles_internal_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW2065_STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_GATE.md"


def main() -> None:
    r2063 = json.loads((ROOT / "report_qw2063_derivational_reconstruction_shared_flavor_basis.json").read_text(encoding="utf-8"))
    r2064 = json.loads((ROOT / "report_qw2064_micro_derived_renormalization_constants_gate.json").read_text(encoding="utf-8"))

    f63: Dict[str, bool] = {k: bool(v) for k, v in r2063["flags"].items()}
    # Override only the foundational strict flag by explicit gate result from QW-2064.
    strict_foundational = bool(str(r2064["verdict"]).startswith("MICRO_DERIVED_RENORMALIZATION_CONSTANTS_GATE_PASS"))

    flags = dict(f63)
    flags["strict_first_principles_foundational_constants_derived"] = strict_foundational

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    physical_keys = {
        "mass_mean_rel_pct_le_max",
        "mass_max_rel_pct_le_max",
        "mass_tau_charm_ratio_err_le_max",
        "ckm_mean_rel_pct_le_max",
        "pmns_mean_rel_pct_le_max",
        "gw_sep_ge_min",
        "gw_adv_ge_min",
        "gw_auc_ge_min",
        "gw_control_gap_le_max",
    }
    physical_pass = bool(all(flags[k] for k in physical_keys))

    strict_internal_pass = bool(all(flags.values()))

    if strict_internal_pass:
        verdict = "STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_PASS"
    elif physical_pass:
        verdict = "STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_PARTIAL"
    else:
        verdict = "STRICT_FIRST_PRINCIPLES_INTERNAL_CLOSURE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "triad_physical": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
            "foundational_gate": "report_qw2064_micro_derived_renormalization_constants_gate.json",
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "physical_pass": physical_pass,
        "strict_internal_pass": strict_internal_pass,
        "foundation_gate_verdict": r2064["verdict"],
        "foundation_ci_warning": bool(r2064.get("ci_warning", False)),
        "verdict": verdict,
        "required_next_step": (
            "RUN_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE"
            if strict_internal_pass
            else "REPAIR_INTERNAL_STRICT_FLAGS"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2065: STRICT FIRST-PRINCIPLES INTERNAL CLOSURE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- physical_pass: {physical_pass}",
        f"- strict_internal_pass: {strict_internal_pass}",
        f"- foundation_gate_verdict: {out['foundation_gate_verdict']}",
        f"- foundation_ci_warning: {out['foundation_ci_warning']}",
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

    print(f"[QW-2065] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2065] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2065] verdict={verdict} pass_count={pass_count}/{total_flags} strict_internal_pass={strict_internal_pass}")


if __name__ == "__main__":
    main()
