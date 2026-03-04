#!/usr/bin/env python3
"""
QW-2104: T3/T4 strict preflight gate.

Purpose:
- merge strict input readiness + provenance pre-gates with downstream T3/T4 gates
  into one deterministic status artifact,
- block accidental overclaim by checking logical consistency.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
R2106 = ROOT / "report_qw2106_strict_external_input_intake_gate.json"
R2099 = ROOT / "report_qw2099_hz_external_decoupling_autocollector.json"
R2102 = ROOT / "report_qw2102_hz_decoupling_identifiability_gate.json"
R2090 = ROOT / "report_qw2090_h0_lambda_decoupling_gate.json"
R2101 = ROOT / "report_qw2101_gnewton_bridge_external_autocollector.json"
R2103 = ROOT / "report_qw2103_gnewton_dimensionless_provenance_gate.json"
R2092 = ROOT / "report_qw2092_gnewton_si_bridge_gate.json"

OUT_JSON = ROOT / "report_qw2104_t3t4_strict_preflight_gate.json"
OUT_MD = ROOT / "RAPORT_QW2104_T3T4_STRICT_PREFLIGHT_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_optional_json(path: Path) -> Dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def verdict_is(v: str, needle: str) -> bool:
    return str(v).strip() == needle


def intake_hz_ready(d2106: Dict | None) -> bool:
    if d2106 is None:
        return False
    f = d2106.get("hz_flags", {})
    needed = [
        "hz_csv_exists",
        "hz_csv_parse_ok",
        "hz_meta_exists",
        "hz_meta_complete",
        "hz_provenance_anchor_free",
        "hz_n_nodes_ge_5",
        "hz_z_span_ge_0p8",
        "hz_e_span_ge_1p0",
        "hz_design_condition_lt_8",
    ]
    return all(bool(f.get(k, False)) for k in needed)


def intake_g_ready(d2106: Dict | None) -> bool:
    if d2106 is None:
        return False
    f = d2106.get("g_flags", {})
    needed = [
        "g_json_exists",
        "g_json_parse_ok",
        "g_meta_exists",
        "g_meta_complete",
        "g_provenance_anchor_free",
        "g_not_seeded_from_registry",
        "g_origin_external_dimensionless",
        "g_dimensionless_present_positive",
        "g_si_not_primary",
    ]
    return all(bool(f.get(k, False)) for k in needed)


def main() -> None:
    d2106 = load_optional_json(R2106)
    d2099 = load_json(R2099)
    d2102 = load_json(R2102)
    d2090 = load_json(R2090)
    d2101 = load_json(R2101)
    d2103 = load_json(R2103)
    d2092 = load_json(R2092)

    intake_hz_strict_ready = intake_hz_ready(d2106)
    intake_g_strict_ready = intake_g_ready(d2106)

    hz_input_ready = bool(d2099.get("strict_ready", False))
    hz_identifiability_pass = verdict_is(
        d2102.get("verdict", ""), "HZ_DECOUPLING_IDENTIFIABILITY_GATE_PASS_STRICT_READY"
    )
    hz_decoupling_strict = verdict_is(d2090.get("verdict", ""), "H0_LAMBDA_DECOUPLING_GATE_PASS_STRICT")

    g_bridge_ready = bool(d2101.get("bridge", {}).get("strict_provenance_ready", False))
    g_provenance_pass = verdict_is(
        d2103.get("verdict", ""), "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY"
    )
    g_si_bridge_strict = verdict_is(d2092.get("verdict", ""), "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT")

    flags = {
        "intake_hz_strict_ready": intake_hz_strict_ready,
        "hz_input_strict_ready": hz_input_ready,
        "hz_identifiability_gate_pass": hz_identifiability_pass,
        "hz_decoupling_gate_strict_pass": hz_decoupling_strict,
        "intake_g_strict_ready": intake_g_strict_ready,
        "g_bridge_input_strict_ready": g_bridge_ready,
        "g_provenance_gate_pass": g_provenance_pass,
        "g_si_bridge_gate_strict_pass": g_si_bridge_strict,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    defects = []
    if not intake_hz_strict_ready and hz_decoupling_strict:
        defects.append("QW-2090 strict pass while QW-2106 H(z) intake strict readiness is false.")
    if not hz_input_ready and hz_decoupling_strict:
        defects.append("QW-2090 strict pass while QW-2099 strict input readiness is false.")
    if not hz_identifiability_pass and hz_decoupling_strict:
        defects.append("QW-2090 strict pass while QW-2102 identifiability gate is not strict-pass.")
    if not intake_g_strict_ready and g_si_bridge_strict:
        defects.append("QW-2092 strict pass while QW-2106 G_newton intake strict readiness is false.")
    if not g_bridge_ready and g_si_bridge_strict:
        defects.append("QW-2092 strict pass while QW-2101 strict provenance readiness is false.")
    if not g_provenance_pass and g_si_bridge_strict:
        defects.append("QW-2092 strict pass while QW-2103 provenance gate is not strict-pass.")

    has_defect = bool(defects)
    all_pass = bool(all(flags.values()))

    if all_pass and not has_defect:
        verdict = "T3T4_STRICT_PREFLIGHT_GATE_PASS"
        required_next_step = "PROMOTE_T3T4_TO_STRICT_UNRESOLVED_REDUCTION_CHAIN"
    elif has_defect:
        verdict = "T3T4_STRICT_PREFLIGHT_GATE_FAIL_LOGIC_DEFECT"
        required_next_step = "FIX_PRECONDITION_LOGIC_CONTRADICTIONS_AND_RERUN"
    else:
        verdict = "T3T4_STRICT_PREFLIGHT_GATE_PENDING"
        required_next_step = (
            "PROVIDE_STRICT_READY_RAW_INPUTS_AND_DIMENSIONLESS_GNEWTON_BRIDGE_THEN_RERUN_QW2106_QW2099_QW2101_QW2102_QW2103_QW2090_QW2092"
        )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2106": R2106.name if d2106 is not None else None,
            "qw2099": R2099.name,
            "qw2102": R2102.name,
            "qw2090": R2090.name,
            "qw2101": R2101.name,
            "qw2103": R2103.name,
            "qw2092": R2092.name,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "defects": defects,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2104: T3T4 STRICT PREFLIGHT GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: `{v}`")
    lines.extend(
        [
            "",
            "## Defects",
            f"- count: `{len(defects)}`",
        ]
    )
    for d in defects:
        lines.append(f"- {d}")
    lines.extend(
        [
            "",
            "## Artifact",
            f"- JSON: `{OUT_JSON.name}`",
            f"- required_next_step: `{required_next_step}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2104] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2104] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2104] verdict={verdict} pass_count={pass_count}/{total_flags} defects={len(defects)}")


if __name__ == "__main__":
    main()
