#!/usr/bin/env python3
"""
QW-2092: G_newton SI-bridge gate (deterministic, no-retune, no-scan).

Goal:
- derive G_newton from a dimensionless gravity bridge observable in a strictly
  auditable way,
- prevent false strict closure when bridge inputs are anchor-seeded/circular.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "gnewton_si_bridge_input_qw2092.json"
OUT_JSON = ROOT / "report_qw2092_gnewton_si_bridge_gate.json"
OUT_MD = ROOT / "RAPORT_QW2092_GNEWTON_SI_BRIDGE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_ref(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    if ref == 0.0:
        return float("inf") if pred != 0.0 else 0.0
    return abs(pred - ref) / abs(ref) * 100.0


def source_metadata_complete(in_data: Dict) -> bool:
    src = str(in_data.get("source", "")).strip()
    citation = str(in_data.get("citation", "")).strip()
    ref_url = str(in_data.get("reference_url", "")).strip()
    src_sha = str(in_data.get("source_sha256", "")).strip()
    src_ver = str(in_data.get("source_version", "")).strip()
    low = f"{src} {citation} {ref_url}".lower()
    not_placeholder = bool(src) and ("placeholder" not in low)
    has_reference = bool(citation or ref_url)
    has_integrity = bool(src_sha or src_ver)
    return bool(not_placeholder and has_reference and has_integrity)


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2092 G_newton SI-bridge gate")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="SI bridge input JSON.")
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    r2073 = load_json(ROOT / "report_qw2073_radiative_channels_closure_upgrade.json")
    groups = reg["groups"]

    g_ref_item = get_ref(groups, "G_newton")
    g_ref = float(g_ref_item["value"])
    g_tol = float(g_ref_item["tolerance_rel_pct"])

    in_path = Path(args.input).resolve()
    input_present = in_path.exists()
    in_data = load_json(in_path) if input_present else {}

    # Default physical constants for SI conversion.
    hbar_si = float(in_data.get("hbar_si", 1.054571817e-34))
    c_si = float(in_data.get("c_si", 299792458.0))
    gev_to_j = float(in_data.get("gev_to_j", 1.602176634e-10))

    candidate_source = "fallback_qw2073_gr_samples_anchor_seeded"
    bridge_failed_reason = ""

    if input_present and in_data.get("g_dimensionless_mu_ref") is not None:
        g_dimless = float(in_data["g_dimensionless_mu_ref"])
        mu_ref_gev = float(in_data.get("mu_ref_gev", 1.0))
        candidate_source = str(in_data.get("source", "input_json"))
        independent_bridge = True
        seeded_from_registry = bool(in_data.get("seeded_from_registry", False))
        bridge_origin = str(in_data.get("bridge_observable_origin", "unknown")).strip()
    else:
        # Deterministic fallback from QW-2073 diagnostics (explicitly anchor-seeded by construction).
        sample = None
        for row in r2073.get("diagnostics", {}).get("gr_samples", []):
            if abs(float(row.get("mu_gev", 0.0)) - 1.0) < 1e-12:
                sample = row
                break
        if sample is None:
            raise RuntimeError("QW-2092 fallback requires mu=1 GeV sample in QW-2073 diagnostics.")
        g_dimless = float(sample["g_dimensionless"])
        mu_ref_gev = 1.0
        independent_bridge = False
        seeded_from_registry = True
        bridge_origin = "qw2073_anchor_seeded_fallback"
        bridge_failed_reason = "independent_bridge_input_missing"

    if g_dimless <= 0.0 or mu_ref_gev <= 0.0:
        raise RuntimeError("Unphysical bridge input: g_dimensionless and mu_ref_gev must be > 0.")

    if bridge_origin == "backsolved_from_g_si" and not bridge_failed_reason:
        bridge_failed_reason = "bridge_observable_backsolved_from_g_si_tautological_for_strict"

    g_nat_gev_m2 = g_dimless / (mu_ref_gev * mu_ref_gev)
    g_si = float(g_nat_gev_m2 * hbar_si * (c_si**5) / (gev_to_j**2))
    g_rel = rel_err_pct(g_si, g_ref)

    flags = {
        "bridge_input_present": bool(input_present),
        "source_metadata_complete": bool(source_metadata_complete(in_data) if input_present else False),
        "independent_bridge_available": bool(independent_bridge),
        "bridge_not_seeded_from_registry": bool(not seeded_from_registry),
        "bridge_not_backsolved_from_g_si": bool(bridge_origin != "backsolved_from_g_si"),
        "conversion_finite_positive": bool(np.isfinite(g_si) and g_si > 0.0),
        "within_registry_tolerance": bool(g_rel <= g_tol),
        "no_anchor_feedback_loop": bool(in_data.get("provenance_anchor_free", False)),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_pass = bool(all(flags.values()))

    if strict_pass:
        verdict = "GNEWTON_SI_BRIDGE_GATE_PASS_STRICT"
        update = {
            "id": "G_newton",
            "predicted_value": float(g_si),
            "status": "derived",
            "strict_level": "strict_internal_gate",
            "method": "qw2092_dimensionless_to_si_bridge_from_independent_observable",
            "notes": "G_newton derived via independent dimensionless bridge observable with explicit SI conversion.",
        }
    else:
        verdict = "GNEWTON_SI_BRIDGE_GATE_PENDING_NONCLOSING"
        update = {
            "id": "G_newton",
            "predicted_value": float(g_si),
            "status": "derived_nofit_anchor_dependent",
            "strict_level": "physical_relation_anchor_dependent",
            "method": "qw2092_fallback_bridge_anchor_seeded",
            "notes": "Bridge remains anchor-dependent/coupled; strict closure blocked until independent bridge observable is provided.",
        }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "radiative_upgrade": "report_qw2073_radiative_channels_closure_upgrade.json",
            "input_json": str(in_path) if input_present else None,
        },
        "bridge": {
            "candidate_source": candidate_source,
            "bridge_failed_reason": bridge_failed_reason,
            "g_dimensionless_mu_ref": g_dimless,
            "mu_ref_gev": mu_ref_gev,
            "bridge_observable_origin": bridge_origin,
            "g_nat_gev_m2": g_nat_gev_m2,
            "g_si": g_si,
        },
        "reference": {
            "g_ref": g_ref,
            "tol_rel_pct": g_tol,
            "rel_err_pct": g_rel,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "update": update,
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_QW2092_IN_QW2069_AND_RERUN_FULL_CLOSURE_CHAIN"
            if strict_pass
            else "PROVIDE_INDEPENDENT_DIMENSIONLESS_BRIDGE_OBSERVABLE_AND_RERUN_QW2092"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2092: GNEWTON SI BRIDGE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- candidate_source: {candidate_source}",
        f"- bridge_failed_reason: {bridge_failed_reason}",
        "",
        "## Derived",
        f"- G_newton candidate: {g_si:.12e} m^3 kg^-1 s^-2",
        f"- rel_err_pct vs registry: {g_rel:.6f} (tol={g_tol:.3f})",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2092] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2092] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2092] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
