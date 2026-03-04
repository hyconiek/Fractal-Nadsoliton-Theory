#!/usr/bin/env python3
"""
QW-2079: PMNS CP external autocollector for QW-2077.

Purpose:
- populate the PMNS observation block in QW-2077 input JSON,
- enforce explicit numeric validation (range and CI ordering),
- preserve provenance metadata.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
DEFAULT_TEMPLATE = ROOT / "empirical_observations_input_qw2077.template.json"
DEFAULT_IN = ROOT / "empirical_observations_input_qw2077.gw_autocollected.json"
DEFAULT_OUT = ROOT / "empirical_observations_input_qw2077.gw_pmns_autocollected.json"
OUT_JSON = ROOT / "report_qw2079_pmns_cp_external_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2079_PMNS_CP_EXTERNAL_AUTOCOLLECTOR.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2079 PMNS CP external autocollector")
    p.add_argument("--input-observation", default=str(DEFAULT_IN), help="Base observation JSON to update.")
    p.add_argument("--output-observation", default=str(DEFAULT_OUT), help="Output observation JSON path.")
    p.add_argument("--sin-delta-central", type=float, required=True, help="Observed central sin(delta_cp_pmns).")
    p.add_argument("--sin-delta-ci95-low", type=float, required=True, help="Observed 95% CI lower bound.")
    p.add_argument("--sin-delta-ci95-high", type=float, required=True, help="Observed 95% CI upper bound.")
    p.add_argument("--source", required=True, help="External source label/citation.")
    args = p.parse_args()

    in_path = Path(args.input_observation).resolve()
    out_path = Path(args.output_observation).resolve()

    if in_path.exists():
        obs = load_json(in_path)
    elif DEFAULT_TEMPLATE.exists():
        obs = load_json(DEFAULT_TEMPLATE)
    else:
        obs = {
            "pmns_cp": {},
            "cosmology_weff_nodes": [],
            "gw_future_holdout": {},
            "notes": "",
        }

    c = float(args.sin_delta_central)
    lo = float(args.sin_delta_ci95_low)
    hi = float(args.sin_delta_ci95_high)
    if not (-1.0 <= c <= 1.0 and -1.0 <= lo <= 1.0 and -1.0 <= hi <= 1.0):
        raise ValueError("sin(delta) values must be within [-1, 1].")
    if lo > hi:
        raise ValueError("CI bounds must satisfy low <= high.")

    obs["pmns_cp"] = {
        "sin_delta_central": c,
        "sin_delta_ci95_low": lo,
        "sin_delta_ci95_high": hi,
        "source": str(args.source),
    }
    note = (
        f"QW-2079 PMNS block updated at {datetime.now(timezone.utc).isoformat()} "
        f"from source={args.source}."
    )
    obs["notes"] = ((obs.get("notes") or "").strip() + " " + note).strip()
    out_path.write_text(json.dumps(obs, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "PMNS_CP_EXTERNAL_AUTOCOLLECTED",
        "input_observation": str(in_path),
        "output_observation": str(out_path),
        "pmns_cp": obs["pmns_cp"],
        "required_next_step": f"RUN_QW2077_WITH:{out_path.name}",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2079: PMNS CP EXTERNAL AUTOCOLLECTOR",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Input observation JSON: `{in_path}`",
        f"- Output observation JSON: `{out_path}`",
        "",
        "## PMNS Inputs",
        f"- sin_delta_central: `{c:.9f}`",
        f"- sin_delta_ci95_low: `{lo:.9f}`",
        f"- sin_delta_ci95_high: `{hi:.9f}`",
        f"- source: `{args.source}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2079] Saved observation JSON: {out_path}")
    print(f"[QW-2079] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2079] Saved report MD:   {OUT_MD.name}")
    print(f"[QW-2079] verdict={out['verdict']}")


if __name__ == "__main__":
    main()

