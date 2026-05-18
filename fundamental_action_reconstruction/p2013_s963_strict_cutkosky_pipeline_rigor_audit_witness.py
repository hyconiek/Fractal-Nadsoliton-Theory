#!/usr/bin/env python3
"""P2013 S963 strict Cutkosky pipeline rigor audit witness.

Audits P1997..P2012 packets for scientific-rigor metadata fields required by the
current continuation discipline: explicit non-closure guard, next-honest-step,
lay explanation, and progress-to-ToE statement.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2013_s963_strict_cutkosky_pipeline_rigor_audit_witness.json"
TS = "2026-05-18T00:00:00+00:00"

TARGETS = [
    "p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json",
    "p1998_s948_strict_cutkosky_channelwise_backend_weight_calibration_witness.json",
    "p1999_s949_strict_cutkosky_backend_kappa_calibrated_channel_solver_witness.json",
    "p2000_s950_strict_cutkosky_loop_kernel_two_channel_witness.json",
    "p2001_s951_strict_cutkosky_full_three_loop_kernel_plus_extra_channel_witness.json",
    "p2002_s952_strict_cutkosky_missing_channel_vs_structural_classifier_witness.json",
    "p2003_s953_strict_cutkosky_classifier_robustness_band_witness.json",
    "p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json",
    "p2005_s955_strict_cutkosky_gx_backend_amplitude_bound_and_classifier_witness.json",
    "p2006_s956_strict_cutkosky_gx_backend_tensor_surrogate_and_covariance_classifier_witness.json",
    "p2007_s957_strict_cutkosky_tensor_vs_scalar_band_comparator_witness.json",
    "p2008_s958_strict_cutkosky_direct_backend_tensor_object_classifier_witness.json",
    "p2009_s959_strict_cutkosky_channelwise_tensor_coupled_covariance_classifier_witness.json",
    "p2011_s961_strict_cutkosky_loop_amplitude_tensor_placeholder_and_coupled_scan_witness.json",
    "p2012_s962_strict_cutkosky_channel_loop_tensor_consistency_audit_witness.json",
]


def load_packet(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def has_text(packet: dict[str, Any], key: str) -> bool:
    v = packet.get(key)
    return isinstance(v, str) and bool(v.strip())


def main() -> None:
    GEN.mkdir(exist_ok=True)

    rows = []
    for name in TARGETS:
        pkt = load_packet(name)
        ok = {
            "has_false_pass_guard": has_text(pkt, "false_pass_guard"),
            "has_next_honest_step": has_text(pkt, "next_honest_step"),
            "has_lay_explanation": has_text(pkt, "lay_explanation"),
            "has_toe_progress_note": has_text(pkt, "toe_progress"),
        }
        rows.append(
            {
                "artifact": name,
                "ledger_id": pkt.get("ledger_id"),
                "checks": ok,
                "all_checks_pass": all(ok.values()),
            }
        )

    passed = sum(1 for r in rows if r["all_checks_pass"])
    total = len(rows)

    out = {
        "ledger_id": "P2013_S963_STRICT_CUTKOSKY_PIPELINE_RIGOR_AUDIT_WITNESS",
        "packet_id": "P2013",
        "stage_id": "S963",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "scope": "P1997..P2012 packet metadata discipline audit",
        "required_fields": [
            "false_pass_guard",
            "next_honest_step",
            "lay_explanation",
            "toe_progress",
        ],
        "audit_table": rows,
        "summary": {"passed": passed, "total": total},
        "gatekeeper_checks": {
            "all_artifacts_present": all(r["ledger_id"] is not None for r in rows),
            "all_required_fields_present": passed == total,
        },
        "result_kind": "PASS_RIGOR_METADATA_AUDIT_WITNESS" if passed == total else "OPEN_RIGOR_METADATA_GAPS",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE" if passed < total else "PASS_WITH_LIMITED_SCOPE",
        "false_pass_guard": "Audit checks metadata discipline only; it does not prove physical closure or strict-core selector discharge.",
        "next_honest_step": "Fill missing toe_progress/discipline fields in failing packets, then rerun P2013 and continue with strict-side non-placeholder loop amplitudes.",
        "toe_progress": "Improves methodological honesty and traceability of strict-lane diagnostics, but does not close ToE or discharge QW-2191.",
        "lay_explanation": "To jest kontrola jakości raportów: sprawdza, czy każdy raport jasno mówi o ograniczeniach, następnym kroku i znaczeniu dla drogi do ToE.",
        "environment": {"python": platform.python_version()},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2013] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
