#!/usr/bin/env python3
"""
QW-2103: G_newton dimensionless provenance gate.

Purpose:
- verify whether the QW-2092 input contains a truly independent dimensionless
  bridge observable (not backsolved from G_SI),
- enforce transparent epistemic labeling before any strict closure claim.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "gnewton_si_bridge_input_qw2092.json"
OUT_JSON = ROOT / "report_qw2103_gnewton_dimensionless_provenance_gate.json"
OUT_MD = ROOT / "RAPORT_QW2103_GNEWTON_DIMENSIONLESS_PROVENANCE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


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
    p = argparse.ArgumentParser(description="QW-2103 G_newton dimensionless provenance gate")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="QW-2092 input JSON")
    args = p.parse_args()

    in_path = Path(args.input).resolve()
    input_present = in_path.exists()
    in_data = load_json(in_path) if input_present else {}

    g_dim = in_data.get("g_dimensionless_mu_ref")
    mu_ref = in_data.get("mu_ref_gev")
    g_si_opt = in_data.get("g_si_input_optional")
    bridge_origin = str(in_data.get("bridge_observable_origin", "")).strip()
    seeded = bool(in_data.get("seeded_from_registry", False))
    anchor_free = bool(in_data.get("provenance_anchor_free", False))

    flags = {
        "input_present": bool(input_present),
        "source_metadata_complete": bool(source_metadata_complete(in_data) if input_present else False),
        "g_dimensionless_present_positive": bool(
            g_dim is not None and np.isfinite(float(g_dim)) and float(g_dim) > 0.0
        ),
        "mu_ref_present_positive": bool(
            mu_ref is not None and np.isfinite(float(mu_ref)) and float(mu_ref) > 0.0
        ),
        "bridge_origin_external_dimensionless": bool(bridge_origin == "external_dimensionless_observable"),
        "not_seeded_from_registry": bool(not seeded),
        "provenance_anchor_free": bool(anchor_free),
        "g_si_not_primary_input": bool(g_si_opt is None),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_ready = bool(all(flags.values()))
    verdict = (
        "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PASS_STRICT_READY"
        if strict_ready
        else "GNEWTON_DIMENSIONLESS_PROVENANCE_GATE_PENDING_NONCLOSING"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {"input_json": str(in_path) if input_present else None},
        "input_summary": {
            "bridge_observable_origin": bridge_origin,
            "seeded_from_registry": seeded,
            "provenance_anchor_free": anchor_free,
            "g_dimensionless_mu_ref": (float(g_dim) if g_dim is not None else None),
            "mu_ref_gev": (float(mu_ref) if mu_ref is not None else None),
            "g_si_input_optional": (float(g_si_opt) if g_si_opt is not None else None),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "RUN_QW2092_WITH_THIS_INPUT_STRICT_READY"
            if strict_ready
            else "PROVIDE_DIRECT_EXTERNAL_DIMENSIONLESS_BRIDGE_OBSERVABLE_AND_RERUN_QW2101_QW2103_QW2092"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2103: GNEWTON DIMENSIONLESS PROVENANCE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Input Summary",
        f"- bridge_observable_origin: `{bridge_origin}`",
        f"- seeded_from_registry: `{seeded}`",
        f"- provenance_anchor_free: `{anchor_free}`",
        f"- g_dimensionless_mu_ref: `{out['input_summary']['g_dimensionless_mu_ref']}`",
        f"- g_si_input_optional: `{out['input_summary']['g_si_input_optional']}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2103] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2103] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2103] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

