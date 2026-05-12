#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _exists(rel: str) -> bool:
    return (FAR / rel).exists()


def main() -> None:
    claims = [
        {
            "id": "C1_gauge_couplings_g_gp_g3",
            "status": "strict_candidate",
            "evidence": ["A6_GAUGE_RECONSTRUCTION_SPEC.md"],
            "residual_artifact": "generated/p1358_residual_table_summary.json",
            "notes": "A6 marks partial strict-core derivation lane; not yet residual-verified match.",
        },
        {
            "id": "C2_weinberg_angle_legacy_formula",
            "status": "legacy_only",
            "evidence": ["P1359_TEX_DERIVATION_REANALYSIS_LEGACY_TO_R81_STRICT_CORRECTION_PACKET_PL.md"],
            "residual_artifact": None,
            "notes": "Legacy identity requires strict successor/provenance theorem before promotion.",
        },
        {
            "id": "C3_fine_structure_successor",
            "status": "strict_candidate",
            "evidence": ["p82_strict_side_role_equivalence_probe_for_legacy_fine_structure_role.py"],
            "residual_artifact": None,
            "notes": "Candidate exists but role-equivalence verdict missing.",
        },
        {
            "id": "C4_mass_proxy_f704_to_gev",
            "status": "nonstrict_proxy",
            "evidence": ["P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"],
            "residual_artifact": None,
            "notes": "Explicitly nonstrict calibration layer.",
        },
        {
            "id": "C5_kernel_only_first_prediction_run",
            "status": "strict_candidate",
            "evidence": ["generated/p1358_kernel_value_generator_summary.json"],
            "residual_artifact": "generated/p1358_residual_table_summary.json",
            "notes": "First non-template run done; current residual verdict FAIL.",
        },
    ]

    for c in claims:
        for e in c["evidence"]:
            c.setdefault("evidence_exists", True)
            c["evidence_exists"] = c["evidence_exists"] and _exists(e)
        ra = c["residual_artifact"]
        c["residual_exists"] = _exists(ra) if ra else False

    summary = {
        "packet": "P1361",
        "as_of": "2026-05-12",
        "claims": claims,
        "status_counts": {
            "strict_verified": sum(1 for c in claims if c["status"] == "strict_verified"),
            "strict_candidate": sum(1 for c in claims if c["status"] == "strict_candidate"),
            "nonstrict_proxy": sum(1 for c in claims if c["status"] == "nonstrict_proxy"),
            "legacy_only": sum(1 for c in claims if c["status"] == "legacy_only"),
        },
        "next_priority": "P1362_RESIDUAL_BENCHMARK_AND_UNCERTAINTY_BUDGET_FOR_STRICT_CANDIDATES",
    }

    out = GEN / "p1361_far_constant_claims_scoreboard_summary.json"
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1361] wrote {out}")


if __name__ == "__main__":
    main()
