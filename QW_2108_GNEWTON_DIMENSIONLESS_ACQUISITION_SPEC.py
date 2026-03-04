#!/usr/bin/env python3
"""
QW-2108: strict acquisition specification for external dimensionless G_newton bridge.

Purpose:
- derive the exact target range for an external dimensionless bridge observable
  (g_dimensionless_mu_ref) from registry G_newton and frozen conversion constants,
- provide a deterministic, auditable acquisition contract for independent teams,
- avoid backsolved/SI-primary circularity in QW-2101/QW-2103/QW-2092 chain.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
R2068 = ROOT / "report_qw2068_sm_gr_parameter_registry.json"
OUT_JSON = ROOT / "report_qw2108_gnewton_dimensionless_acquisition_spec.json"
OUT_MD = ROOT / "RAPORT_QW2108_GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_registry_param(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2108 strict external G bridge acquisition spec")
    p.add_argument("--mu-ref-gev", type=float, default=1.0, help="Reference scale in GeV for bridge observable")
    p.add_argument("--hbar-si", type=float, default=1.054571817e-34, help="Reduced Planck constant (SI)")
    p.add_argument("--c-si", type=float, default=299792458.0, help="Speed of light (SI)")
    p.add_argument("--gev-to-j", type=float, default=1.602176634e-10, help="GeV to Joule conversion")
    args = p.parse_args()

    mu_ref = float(args.mu_ref_gev)
    hbar_si = float(args.hbar_si)
    c_si = float(args.c_si)
    gev_to_j = float(args.gev_to_j)

    if not (math.isfinite(mu_ref) and mu_ref > 0.0):
        raise ValueError("mu_ref_gev must be finite and > 0")
    if not all(math.isfinite(x) and x > 0.0 for x in [hbar_si, c_si, gev_to_j]):
        raise ValueError("hbar_si/c_si/gev_to_j must be finite and > 0")

    reg = load_json(R2068)
    g_item = get_registry_param(reg["groups"], "G_newton")
    g_ref = float(g_item["value"])
    tol_pct = float(g_item["tolerance_rel_pct"])
    tol_rel = tol_pct / 100.0

    # Inverse of QW-2101 bridge conversion.
    g_dim_target = g_ref * (mu_ref * mu_ref) * (gev_to_j**2) / (hbar_si * (c_si**5))
    g_dim_min = g_dim_target * (1.0 - tol_rel)
    g_dim_max = g_dim_target * (1.0 + tol_rel)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": R2068.name,
            "registry_parameter_id": "G_newton",
        },
        "registry_target": {
            "parameter_id": "G_newton",
            "g_ref_si": g_ref,
            "g_ref_unit": "m^3 kg^-1 s^-2",
            "tolerance_rel_pct": tol_pct,
        },
        "bridge_spec": {
            "mu_ref_gev": mu_ref,
            "conversion_constants": {
                "hbar_si": hbar_si,
                "c_si": c_si,
                "gev_to_j": gev_to_j,
            },
            "g_dimensionless_target": g_dim_target,
            "g_dimensionless_acceptance_range": {
                "min": g_dim_min,
                "max": g_dim_max,
            },
        },
        "strict_contract": {
            "required_bridge_observable_origin": "external_dimensionless_observable",
            "must_not_be_backsolved_from_g_si": True,
            "must_not_be_seeded_from_registry": True,
            "provenance_anchor_free_required": True,
            "g_si_input_optional_must_be_null_for_strict": True,
            "required_metadata_fields": [
                "source_dataset_name",
                "citation",
                "reference_url",
                "source_version",
                "provenance_anchor_free",
                "seeded_from_registry",
                "bridge_observable_origin",
            ],
        },
        "output_contract_examples": {
            "external_gnewton_bridge_qw2101.json": {
                "mu_ref_gev": mu_ref,
                "g_dimensionless_mu_ref": "MEASURED_EXTERNAL_VALUE_IN_RANGE",
                "bridge_observable_origin": "external_dimensionless_observable",
                "notes": "Do not backsolve from G_SI; provide direct external observable pipeline output."
            },
            "external_gnewton_bridge_qw2101.metadata.json": {
                "source_dataset_name": "REAL_EXTERNAL_DATASET_NAME",
                "citation": "REAL_CITATION",
                "reference_url": "REAL_URL_OR_DOI",
                "source_version": "REAL_RELEASE_TAG",
                "provenance_anchor_free": True,
                "seeded_from_registry": False,
                "bridge_observable_origin": "external_dimensionless_observable",
            },
        },
        "verdict": "GNEWTON_DIMENSIONLESS_ACQUISITION_SPEC_READY",
        "required_next_step": "COLLECT_EXTERNAL_DIMENSIONLESS_BRIDGE_OBSERVABLE_WITH_STRICT_METADATA_THEN_RERUN_QW2106_QW2101_QW2103_QW2092",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2108: GNEWTON DIMENSIONLESS ACQUISITION SPEC",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Registry Target",
        f"- G_ref (SI): `{g_ref:.12e}`",
        f"- tolerance_rel_pct: `{tol_pct:.6g}`",
        "",
        "## Dimensionless Bridge Spec",
        f"- mu_ref_gev: `{mu_ref:.12g}`",
        f"- g_dimensionless_target: `{g_dim_target:.12e}`",
        f"- acceptance_range: `[{g_dim_min:.12e}, {g_dim_max:.12e}]`",
        "",
        "## Strict Contract",
        "- bridge_observable_origin must be `external_dimensionless_observable`",
        "- value must NOT be backsolved from G_SI",
        "- provenance_anchor_free must be true",
        "- seeded_from_registry must be false",
        "- g_si_input_optional must be null for strict path",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2108] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2108] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2108] verdict={out['verdict']} g_dim_target={g_dim_target:.12e}")


if __name__ == "__main__":
    main()
