#!/usr/bin/env python3
"""
QW-2015: True external beta-channel v2 readiness + information-richness gate.

Validates package assembled by QW-2014 and certifies whether it is suitable
for strict blind intervention evaluation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2015_true_external_beta_channel_v2_readiness_gate.json"
OUT_MD = ROOT / "RAPORT_QW2015_TRUE_EXTERNAL_BETA_CHANNEL_V2_READINESS_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    base = ROOT / "external_confirmatory_v2" / "beta_channel_true_external_v2"
    pairs = base / "beta_channel_pairs.csv"
    events = base / "intervention_events.csv"
    freeze = base / "protocol_freeze.json"
    manifest = base / "manifest_beta_channel.json"

    missing_files: List[str] = []
    for p in [pairs, events, freeze, manifest]:
        if not p.exists():
            missing_files.append(str(p))

    if missing_files:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_dir": str(base),
            "missing_files": missing_files,
            "verdict": "TRUE_EXTERNAL_BETA_CHANNEL_V2_NOT_READY_MISSING_FILES",
            "readiness": "NOT_READY",
            "required_next_step": "RERUN_QW2014_AUTOCOLLECTOR_V2",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "\n".join(
                [
                    "# RAPORT QW-2015: TRUE EXTERNAL BETA-CHANNEL V2 READINESS GATE",
                    "",
                    f"- Data UTC: {out['generated_utc']}",
                    f"- Readiness: **{out['readiness']}**",
                    f"- Verdict: **{out['verdict']}**",
                    "",
                    "## Missing Files",
                    *[f"- {x}" for x in missing_files],
                    "",
                    "## Required Next Step",
                    f"- {out['required_next_step']}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        print(f"[QW-2015] verdict={out['verdict']}")
        return

    d_manifest = load_json(manifest)
    d_freeze = load_json(freeze)
    d2014 = load_json(ROOT / "report_qw2014_true_external_beta_channel_autocollector_v2.json")

    # Manifest checks.
    provider = str(d_manifest.get("dataset", {}).get("provider", ""))
    ext_stmt = str(d_manifest.get("dataset", {}).get("externality_statement", "")).lower()
    files = d_manifest.get("files", [])
    roles = {str(x.get("role", "")) for x in files if isinstance(x, dict)}

    # Data schema checks.
    df_pairs = pd.read_csv(pairs)
    df_events = pd.read_csv(events)

    req_pairs = [
        "pair_id",
        "theta_deg",
        "hxy",
        "f_std",
        "f_autoc1",
        "f_switch",
        "f_slope",
        "intervention_id",
        "regime",
    ]
    req_events = [
        "intervention_id",
        "intervention_type",
        "source_reference",
        "start_utc",
        "end_utc",
        "is_preregistered",
        "is_exogenous",
    ]

    missing_pair_cols = [c for c in req_pairs if c not in df_pairs.columns]
    missing_event_cols = [c for c in req_events if c not in df_events.columns]

    n_pairs_total = int(len(df_pairs))
    n_events_total = int(len(df_events))
    n_pre = int((df_pairs["regime"].astype(str).str.lower() == "pre").sum())
    n_post = int((df_pairs["regime"].astype(str).str.lower() == "post").sum())

    # Information-richness diagnostics.
    frac_f_std_eq0 = float((df_pairs["f_std"].to_numpy() == 0.0).mean())
    frac_f_autoc1_eq0 = float((df_pairs["f_autoc1"].to_numpy() == 0.0).mean())
    frac_f_switch_eq0 = float((df_pairs["f_switch"].to_numpy() == 0.0).mean())
    frac_f_slope_eq0 = float((df_pairs["f_slope"].to_numpy() == 0.0).mean())

    interventions_used = int(df_pairs["intervention_id"].astype(str).nunique())
    median_abs_autoc1 = float(df_pairs["f_autoc1"].abs().median())
    median_std = float(df_pairs["f_std"].median())

    neigh_stats = d2014.get("feature_richness", {})
    median_neigh = float(neigh_stats.get("median_local_neigh_n", 0.0))

    flags = {
        "externality_ok": bool(
            provider == "NANOGRAV_GWOSC_EXTERNAL_PUBLIC"
            and ("independent" in ext_stmt)
            and ("external" in ext_stmt)
            and ("internal_proxy" not in ext_stmt)
        ),
        "manifest_roles_ok": bool({"beta_pairs", "intervention_events", "protocol_freeze"}.issubset(roles)),
        "schema_pairs_ok": bool(len(missing_pair_cols) == 0),
        "schema_events_ok": bool(len(missing_event_cols) == 0),
        "regime_split_ok": bool(n_pre > 0 and n_post > 0),
        "power_floor_ok": bool(n_pairs_total >= 1200),
        "feature_non_degenerate_ok": bool(
            frac_f_std_eq0 <= 0.10
            and frac_f_slope_eq0 <= 0.10
            and median_neigh >= 5.0
            and median_abs_autoc1 >= 0.05
            and median_std >= 0.01
        ),
        "intervention_diversity_ok": bool(interventions_used >= 3),
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "TRUE_EXTERNAL_BETA_CHANNEL_V2_READY_STRICT"
        readiness = "READY"
        required_next_step = "RUN_V2_BLIND_INTERVENTION_AND_TRIAD_VALIDATION"
    else:
        verdict = "TRUE_EXTERNAL_BETA_CHANNEL_V2_NOT_READY"
        readiness = "NOT_READY"
        required_next_step = "REPAIR_V2_PACKAGE_AND_RERUN_GATE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_dir": str(base),
        "counts": {
            "n_pairs_total": n_pairs_total,
            "n_events_total": n_events_total,
            "n_pre": n_pre,
            "n_post": n_post,
            "interventions_used": interventions_used,
        },
        "missing_pair_columns": missing_pair_cols,
        "missing_event_columns": missing_event_cols,
        "feature_diagnostics": {
            "frac_f_std_eq0": frac_f_std_eq0,
            "frac_f_autoc1_eq0": frac_f_autoc1_eq0,
            "frac_f_switch_eq0": frac_f_switch_eq0,
            "frac_f_slope_eq0": frac_f_slope_eq0,
            "median_abs_autoc1": median_abs_autoc1,
            "median_f_std": median_std,
            "median_local_neigh_n": median_neigh,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2015: TRUE EXTERNAL BETA-CHANNEL V2 READINESS GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Counts",
        f"- n_pairs_total: {n_pairs_total}",
        f"- n_events_total: {n_events_total}",
        f"- n_pre/n_post: {n_pre}/{n_post}",
        f"- interventions_used: {interventions_used}",
        "",
        "## Feature Diagnostics",
        f"- frac_f_std_eq0: {frac_f_std_eq0:.4f}",
        f"- frac_f_autoc1_eq0: {frac_f_autoc1_eq0:.4f}",
        f"- frac_f_switch_eq0: {frac_f_switch_eq0:.4f}",
        f"- frac_f_slope_eq0: {frac_f_slope_eq0:.4f}",
        f"- median_abs_autoc1: {median_abs_autoc1:.4f}",
        f"- median_f_std: {median_std:.4f}",
        f"- median_local_neigh_n: {median_neigh:.2f}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next_step}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2015] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2015] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2015] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
