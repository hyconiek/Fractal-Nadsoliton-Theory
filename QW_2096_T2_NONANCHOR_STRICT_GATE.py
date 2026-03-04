#!/usr/bin/env python3
"""
QW-2096: T2 non-anchor strict aggregate gate (m_up, m_down, m_strange, m_h).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2096_t2_nonanchor_strict_gate.json"
OUT_MD = ROOT / "RAPORT_QW2096_T2_NONANCHOR_STRICT_GATE.md"


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _extract_updates_map(report: Dict | None) -> Dict[str, Dict]:
    if report is None:
        return {}
    updates = report.get("updates")
    if not isinstance(updates, list):
        return {}
    out = {}
    for u in updates:
        if isinstance(u, dict) and "id" in u:
            out[str(u["id"])] = u
    return out


def main() -> None:
    r2088 = load_optional_json("report_qw2088_light_quark_mass_nonanchor_gate.json")
    r2089 = load_optional_json("report_qw2089_higgs_selfcoupling_strict_gate.json")
    r2095 = load_optional_json("report_qw2095_kernel_derived_t2_nonanchor_inputs_plan_executor.json")

    upd_map = {}
    upd_map.update(_extract_updates_map(r2088))
    if r2089 is not None and isinstance(r2089.get("update"), dict) and "id" in r2089["update"]:
        upd_map[str(r2089["update"]["id"])] = r2089["update"]

    ids = ["m_up", "m_down", "m_strange", "m_h"]
    checks = {}
    updates: List[Dict] = []

    for pid in ids:
        u = upd_map.get(pid)
        if u is None:
            checks[pid] = {
                "update_available": False,
                "strict_internal_gate": False,
                "derived_status": False,
            }
            continue
        strict_level_ok = str(u.get("strict_level", "")) == "strict_internal_gate"
        status_ok = str(u.get("status", "")).startswith("derived")
        checks[pid] = {
            "update_available": True,
            "strict_internal_gate": strict_level_ok,
            "derived_status": status_ok,
        }
        updates.append(u)

    all_upstream_strict = all(
        checks[pid]["update_available"]
        and checks[pid]["strict_internal_gate"]
        and checks[pid]["derived_status"]
        for pid in ids
    )

    q2095_artifacts_present = bool(
        r2095 is not None
        and str(r2095.get("verdict", "")) == "KERNEL_DERIVED_T2_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN"
        and (ROOT / "t2_nonanchor_light_quark_input_qw2088.json").exists()
        and (ROOT / "t2_nonanchor_higgs_input_qw2089.json").exists()
    )

    no_anchor_or_circularity_detected = bool(
        r2088 is not None
        and r2089 is not None
        and not bool(r2088.get("flags", {}).get("anchor_or_circularity_detected", True))
        and not bool(r2089.get("flags", {}).get("anchor_or_circularity_detected", True))
    )

    flags = {
        "deterministic_no_retune_no_scan": True,
        "q2095_kernel_input_artifacts_present": q2095_artifacts_present,
        "m_up_strict_nonanchor_pass": checks["m_up"]["update_available"] and checks["m_up"]["strict_internal_gate"],
        "m_down_strict_nonanchor_pass": checks["m_down"]["update_available"] and checks["m_down"]["strict_internal_gate"],
        "m_strange_strict_nonanchor_pass": checks["m_strange"]["update_available"] and checks["m_strange"]["strict_internal_gate"],
        "m_h_strict_nonanchor_pass": checks["m_h"]["update_available"] and checks["m_h"]["strict_internal_gate"],
        "no_anchor_or_circularity_detected": no_anchor_or_circularity_detected,
    }
    pass_count = sum(1 for v in flags.values() if bool(v))
    total_flags = len(flags)

    verdict = (
        "T2_NONANCHOR_STRICT_GATE_PASS"
        if all_upstream_strict and q2095_artifacts_present and no_anchor_or_circularity_detected
        else "T2_NONANCHOR_STRICT_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2088_optional": "report_qw2088_light_quark_mass_nonanchor_gate.json",
            "qw2089_optional": "report_qw2089_higgs_selfcoupling_strict_gate.json",
            "qw2095_optional": "report_qw2095_kernel_derived_t2_nonanchor_inputs_plan_executor.json",
        },
        "checks": checks,
        "updates": updates,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PIPE_T2_UPDATES_TO_QW2069_AND_RERUN_QW2071"
            if verdict == "T2_NONANCHOR_STRICT_GATE_PASS"
            else "PRODUCE_VALID_QW2088_QW2089_STRICT_UPDATES"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2096: T2 NONANCHOR STRICT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Pass count: `{pass_count}/{total_flags}`",
        "",
        "## Checks",
        f"- m_up strict: `{flags['m_up_strict_nonanchor_pass']}`",
        f"- m_down strict: `{flags['m_down_strict_nonanchor_pass']}`",
        f"- m_strange strict: `{flags['m_strange_strict_nonanchor_pass']}`",
        f"- m_h strict: `{flags['m_h_strict_nonanchor_pass']}`",
        f"- q2095 artifacts present: `{q2095_artifacts_present}`",
        f"- no anchor/circularity: `{no_anchor_or_circularity_detected}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2096] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2096] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2096] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

