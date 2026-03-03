#!/usr/bin/env python3
"""
QW-1927: True external beta-channel readiness gate.

Checks whether the required true external intervention package is present and valid
against QW-1925 collection specification.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1927_true_external_beta_channel_readiness_gate.json"
OUT_MD = ROOT / "RAPORT_QW1927_TRUE_EXTERNAL_BETA_CHANNEL_READINESS_GATE.md"


def externality_ok(provider: str, ext_stmt: str, forbidden: List[str], required_tokens: List[str]) -> bool:
    p = provider.upper()
    e = ext_stmt.lower()
    if any(tok.upper() in p for tok in forbidden):
        return False
    if any(tok.lower() in e for tok in ["internal proxy", "not independent"]):
        return False
    return all(tok.lower() in e for tok in required_tokens)


def main() -> None:
    d1925 = json.loads((ROOT / "report_qw1925_true_external_beta_channel_collection_spec.json").read_text(encoding="utf-8"))
    spec = d1925["spec"]

    base = ROOT / "external_confirmatory_v2" / "beta_channel_true_external"

    required_files = list(spec["hard_requirements"]["required_files"])
    missing = [f for f in required_files if not (base / f).exists()]

    if missing:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_dir": str(base),
            "missing_files": missing,
            "verdict": "TRUE_EXTERNAL_BETA_CHANNEL_BLOCKED_MISSING_PACKAGE",
            "readiness": "NOT_READY",
            "required_next_step": "PROVIDE_COMPLETE_TRUE_EXTERNAL_BETA_CHANNEL_PACKAGE",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

        lines = [
            "# RAPORT QW-1927: TRUE EXTERNAL BETA-CHANNEL READINESS GATE",
            "",
            f"- Data UTC: {out['generated_utc']}",
            f"- Verdict: **{out['verdict']}**",
            f"- Readiness: **{out['readiness']}**",
            "",
            "## Missing Files",
        ]
        for f in missing:
            lines.append(f"- {f}")
        lines.extend(
            [
                "",
                "## Required Next Step",
                f"- {out['required_next_step']}",
                "",
                "## Artifacts",
                f"- JSON: `{OUT_JSON.name}`",
            ]
        )
        OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

        print(f"[QW-1927] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1927] Saved MD:   {OUT_MD.name}")
        print(f"[QW-1927] verdict={out['verdict']}")
        return

    manifest = json.loads((base / "manifest_beta_channel.json").read_text(encoding="utf-8"))
    provider = str(manifest.get("dataset", {}).get("provider", ""))
    ext_stmt = str(manifest.get("dataset", {}).get("externality_statement", ""))

    forbidden = list(spec["hard_requirements"]["forbidden_provider_tokens"])
    required_tokens = list(spec["hard_requirements"]["required_externality_statement_tokens"])
    ext_ok = externality_ok(provider, ext_stmt, forbidden=forbidden, required_tokens=required_tokens)

    role_map = {x.get("role"): x for x in manifest.get("files", []) if isinstance(x, dict)}
    role_required = list(spec["hard_requirements"]["required_roles_in_manifest"])
    missing_roles = [r for r in role_required if r not in role_map]

    if missing_roles:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_dir": str(base),
            "externality_ok": ext_ok,
            "missing_manifest_roles": missing_roles,
            "verdict": "TRUE_EXTERNAL_BETA_CHANNEL_BLOCKED_MANIFEST_ROLES",
            "readiness": "NOT_READY",
            "required_next_step": "FIX_MANIFEST_ROLES_AND_RETRY",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

        lines = [
            "# RAPORT QW-1927: TRUE EXTERNAL BETA-CHANNEL READINESS GATE",
            "",
            f"- Data UTC: {out['generated_utc']}",
            f"- Verdict: **{out['verdict']}**",
            f"- Readiness: **{out['readiness']}**",
            f"- externality_ok: {ext_ok}",
            "",
            "## Missing Manifest Roles",
        ]
        for r in missing_roles:
            lines.append(f"- {r}")
        lines.extend(
            [
                "",
                "## Required Next Step",
                f"- {out['required_next_step']}",
                "",
                "## Artifacts",
                f"- JSON: `{OUT_JSON.name}`",
            ]
        )
        OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

        print(f"[QW-1927] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1927] Saved MD:   {OUT_MD.name}")
        print(f"[QW-1927] verdict={out['verdict']}")
        return

    beta_path = base / str(role_map["beta_pairs"]["path"])
    events_path = base / str(role_map["intervention_events"]["path"])

    if not beta_path.exists() or not events_path.exists():
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "candidate_dir": str(base),
            "externality_ok": ext_ok,
            "beta_pairs_exists": bool(beta_path.exists()),
            "intervention_events_exists": bool(events_path.exists()),
            "verdict": "TRUE_EXTERNAL_BETA_CHANNEL_BLOCKED_FILE_LINKS",
            "readiness": "NOT_READY",
            "required_next_step": "FIX_MANIFEST_FILE_LINKS_AND_RETRY",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1927: TRUE EXTERNAL BETA-CHANNEL READINESS GATE\n\n"
            f"- Data UTC: {out['generated_utc']}\n"
            f"- Verdict: **{out['verdict']}**\n"
            f"- Readiness: **{out['readiness']}**\n"
            f"- externality_ok: {ext_ok}\n"
            f"- beta_pairs_exists: {out['beta_pairs_exists']}\n"
            f"- intervention_events_exists: {out['intervention_events_exists']}\n"
            f"\n## Required Next Step\n- {out['required_next_step']}\n",
            encoding="utf-8",
        )
        print(f"[QW-1927] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1927] Saved MD:   {OUT_MD.name}")
        print(f"[QW-1927] verdict={out['verdict']}")
        return

    beta_df = pd.read_csv(beta_path)
    events_df = pd.read_csv(events_path)

    req_pairs = list(spec["dataset_schema"]["beta_channel_pairs_required_columns"])
    req_events = list(spec["dataset_schema"]["intervention_events_required_columns"])

    miss_pairs = [c for c in req_pairs if c not in beta_df.columns]
    miss_events = [c for c in req_events if c not in events_df.columns]

    regime_ok = set(beta_df.get("regime", pd.Series(dtype=str)).astype(str).str.lower().unique()) <= set(
        [x.lower() for x in spec["dataset_schema"]["regime_allowed_values"]]
    )

    pre_count = int((beta_df["regime"].astype(str).str.lower() == "pre").sum()) if "regime" in beta_df.columns else 0
    post_count = int((beta_df["regime"].astype(str).str.lower() == "post").sum()) if "regime" in beta_df.columns else 0

    n_holdout_req = int(spec["power_targets"]["n_holdout_min_for_power_0p90"])
    n_pairs_total = int(len(beta_df))

    hard_flags = {
        "externality_ok": bool(ext_ok),
        "schema_pairs_ok": bool(len(miss_pairs) == 0),
        "schema_events_ok": bool(len(miss_events) == 0),
        "regime_values_ok": bool(regime_ok),
        "has_pre_and_post": bool(pre_count > 0 and post_count > 0),
        "n_pairs_ge_power_target": bool(n_pairs_total >= 2 * n_holdout_req),
    }

    if all(hard_flags.values()):
        verdict = "TRUE_EXTERNAL_BETA_CHANNEL_READY"
        readiness = "READY"
        required_next = "RUN_FROZEN_BLIND_INTERVENTION_EVALUATION"
    else:
        verdict = "TRUE_EXTERNAL_BETA_CHANNEL_NOT_READY"
        readiness = "NOT_READY"
        required_next = "FIX_FAILED_FLAGS_AND_RETRY_READINESS_GATE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_dir": str(base),
        "n_pairs_total": n_pairs_total,
        "n_events_total": int(len(events_df)),
        "n_pre": pre_count,
        "n_post": post_count,
        "power_target_n_holdout": n_holdout_req,
        "missing_pair_columns": miss_pairs,
        "missing_event_columns": miss_events,
        "flags": hard_flags,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1927: TRUE EXTERNAL BETA-CHANNEL READINESS GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- n_pairs_total: {n_pairs_total}",
        f"- n_events_total: {len(events_df)}",
        f"- n_pre/n_post: {pre_count}/{post_count}",
        "",
        "## Flags",
    ]
    for k, v in hard_flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1927] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1927] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1927] verdict={verdict}")


if __name__ == "__main__":
    main()
