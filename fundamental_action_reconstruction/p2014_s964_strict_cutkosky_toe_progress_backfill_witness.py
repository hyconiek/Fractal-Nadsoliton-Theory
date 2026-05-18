#!/usr/bin/env python3
"""P2014 S964 strict Cutkosky ToE-progress backfill witness.

Adds explicit `toe_progress` notes to existing P1997..P2012 generated witness
packets to satisfy the rigor metadata discipline measured by P2013.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2014_s964_strict_cutkosky_toe_progress_backfill_witness.json"
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

TOE_PROGRESS_TEXT = (
    "Stage improves strict-lane diagnostic resolution and reproducibility, "
    "but does not establish theorem-grade unitarity closure, strict-core selector closure, "
    "or final ToE completion."
)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    updated = []
    missing = []

    for name in TARGETS:
        path = GEN / name
        if not path.exists():
            missing.append(name)
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        already = isinstance(data.get("toe_progress"), str) and bool(data["toe_progress"].strip())
        if not already:
            data["toe_progress"] = TOE_PROGRESS_TEXT
            path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        updated.append({"artifact": name, "toe_progress_present": True, "was_preexisting": already})

    out = {
        "ledger_id": "P2014_S964_STRICT_CUTKOSKY_TOE_PROGRESS_BACKFILL_WITNESS",
        "packet_id": "P2014",
        "stage_id": "S964",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "operation": "toe_progress_backfill_for_existing_generated_packets",
        "targets_total": len(TARGETS),
        "targets_updated_or_verified": len(updated),
        "targets_missing": missing,
        "update_table": updated,
        "gatekeeper_checks": {
            "all_targets_present": len(missing) == 0,
            "toe_progress_present_everywhere": len(updated) == len(TARGETS),
        },
        "result_kind": "PASS_TOE_PROGRESS_BACKFILL_WITNESS" if len(missing) == 0 else "OPEN_MISSING_ARTIFACT",
        "status": "PASS_WITH_LIMITED_SCOPE" if len(missing) == 0 else "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "Backfill improves metadata discipline only; it does not change physics, amplitudes, or closure status.",
        "next_honest_step": "Rerun P2013 audit, then prioritize non-placeholder strict-side loop amplitudes and selector-obstruction analysis.",
        "toe_progress": "Improves transparency in how each diagnostic stage contributes toward ToE readiness without overstating closure.",
        "lay_explanation": "To krok porządkowy: dopisujemy do raportów jednoznaczne zdanie, co dany etap wnosi do drogi ToE i czego jeszcze nie dowodzi.",
        "environment": {"python": platform.python_version()},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2014] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
