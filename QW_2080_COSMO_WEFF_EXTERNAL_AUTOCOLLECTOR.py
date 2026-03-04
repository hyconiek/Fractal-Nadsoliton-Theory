#!/usr/bin/env python3
"""
QW-2080: Cosmology w_eff(z) external autocollector for QW-2077.

Purpose:
- populate the cosmology_weff_nodes block in QW-2077 input JSON,
- enforce preregistered node coverage and numeric checks,
- keep provenance metadata for external source tracking.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd


ROOT = Path(__file__).resolve().parent
PREREG_JSON = ROOT / "report_qw2076_empirical_prediction_preregistration.json"
DEFAULT_TEMPLATE = ROOT / "empirical_observations_input_qw2077.template.json"
DEFAULT_IN = ROOT / "empirical_observations_input_qw2077.gw_pmns_autocollected.json"
DEFAULT_OUT = ROOT / "empirical_observations_input_qw2077.full_autocollected.json"
OUT_JSON = ROOT / "report_qw2080_cosmo_weff_external_autocollector.json"
OUT_MD = ROOT / "RAPORT_QW2080_COSMO_WEFF_EXTERNAL_AUTOCOLLECTOR.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def target_z_nodes_from_prereg() -> List[float]:
    data = load_json(PREREG_JSON)
    for pred in data.get("predictions", []):
        if pred.get("id") == "PRED-2076-COSMO-WEFF":
            return [float(x["z"]) for x in pred.get("predicted_nodes", [])]
    raise KeyError("Missing PRED-2076-COSMO-WEFF in preregistration.")


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2080 cosmology w_eff external autocollector")
    p.add_argument("--nodes-csv", required=True, help="CSV with columns: z,w_eff_central,sigma_total[,source].")
    p.add_argument("--input-observation", default=str(DEFAULT_IN), help="Base observation JSON to update.")
    p.add_argument("--output-observation", default=str(DEFAULT_OUT), help="Output observation JSON path.")
    p.add_argument("--source", default="", help="Fallback source if CSV has no 'source' column.")
    args = p.parse_args()

    csv_path = Path(args.nodes_csv).resolve()
    in_path = Path(args.input_observation).resolve()
    out_path = Path(args.output_observation).resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"nodes CSV not found: {csv_path}")

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

    df = pd.read_csv(csv_path)
    need = ["z", "w_eff_central", "sigma_total"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"CSV missing required columns: {miss}")

    z_nodes = target_z_nodes_from_prereg()
    rows: List[Dict] = []
    for z in z_nodes:
        sub = df[(df["z"].astype(float) - z).abs() < 1e-9]
        if len(sub) != 1:
            raise ValueError(f"Expected exactly 1 row for z={z}, found {len(sub)}.")
        r = sub.iloc[0]
        w = float(r["w_eff_central"])
        s = float(r["sigma_total"])
        if s <= 0.0:
            raise ValueError(f"sigma_total must be > 0 for z={z}.")
        src = str(r["source"]) if "source" in df.columns and pd.notna(r["source"]) else str(args.source)
        rows.append(
            {
                "z": float(z),
                "w_eff_central": w,
                "sigma_total": s,
                "source": src,
            }
        )

    obs["cosmology_weff_nodes"] = rows
    note = (
        f"QW-2080 cosmology block updated at {datetime.now(timezone.utc).isoformat()} "
        f"from nodes_csv={csv_path.name}."
    )
    obs["notes"] = ((obs.get("notes") or "").strip() + " " + note).strip()
    out_path.write_text(json.dumps(obs, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "COSMO_WEFF_EXTERNAL_AUTOCOLLECTED",
        "nodes_csv": str(csv_path),
        "input_observation": str(in_path),
        "output_observation": str(out_path),
        "z_nodes": z_nodes,
        "n_nodes": len(rows),
        "required_next_step": f"RUN_QW2077_WITH:{out_path.name}",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2080: COSMO W_EFF EXTERNAL AUTOCOLLECTOR",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Nodes CSV: `{csv_path}`",
        f"- Input observation JSON: `{in_path}`",
        f"- Output observation JSON: `{out_path}`",
        f"- n_nodes: {len(rows)}",
        "",
        "## Nodes",
    ]
    for r in rows:
        lines.append(
            f"- z={r['z']:.6g}, w_eff={r['w_eff_central']:.9f}, sigma={r['sigma_total']:.9f}, source={r['source']}"
        )
    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2080] Saved observation JSON: {out_path}")
    print(f"[QW-2080] Saved report JSON: {OUT_JSON.name}")
    print(f"[QW-2080] Saved report MD:   {OUT_MD.name}")
    print(f"[QW-2080] verdict={out['verdict']} n_nodes={len(rows)}")


if __name__ == "__main__":
    main()

