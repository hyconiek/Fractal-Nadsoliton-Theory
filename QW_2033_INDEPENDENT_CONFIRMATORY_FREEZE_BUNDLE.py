#!/usr/bin/env python3
"""
QW-2033: Build frozen manifest bundle for independent confirmatory replication.

No legacy files are modified. The script creates a new bundle directory with:
- immutable file manifest (sha256 + size),
- minimal runbook for independent rerun of the combined branch.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
BUNDLE_DIR = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2033"
OUT_JSON = ROOT / "report_qw2033_independent_confirmatory_freeze_bundle.json"
OUT_MD = ROOT / "RAPORT_QW2033_INDEPENDENT_CONFIRMATORY_FREEZE_BUNDLE.md"
DATA_SOURCES_DOC = "DATA_SOURCES_EXTERNAL_DOWNLOADS.md"

EXTERNAL_REQUIRED_SOURCES = [
    {
        "name": "NANOGrav 15yr timing archive",
        "dataset": "NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
        "url": "https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1",
        "local_example_path": "external_data/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
        "note": "Keep outside git; pass via --nanograv-archive in autocollector scripts.",
    },
    {
        "name": "GWOSC GWTC event catalog API",
        "dataset": "GWTC event metadata (JSON)",
        "url": "https://www.gw-openscience.org/eventapi/json/GWTC/",
        "local_example_path": "external_data/gwosc_gwtc_eventapi.json",
        "note": "Used as external intervention-event source for beta-channel builds.",
    },
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_file_entry(path: Path) -> Dict:
    if not path.exists():
        return {
            "path": str(path.relative_to(ROOT)),
            "exists": False,
            "size_bytes": None,
            "sha256": None,
        }
    return {
        "path": str(path.relative_to(ROOT)),
        "exists": True,
        "size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
    }


def write_runbook(bundle_dir: Path, files: List[Dict], external_sources: List[Dict]) -> Path:
    runbook = bundle_dir / "RUNBOOK_QW2033.md"
    lines = [
        "# RUNBOOK QW-2033: Independent Confirmatory Replication",
        "",
        "## Scope",
        "- Reproduce combined-branch closure and external preconfirmatory status from frozen artifacts.",
        "- Do not retune sector-specific parameters between mass/flavor/GW.",
        "",
        "## Required Checks",
        "1. Verify SHA256 hashes against `manifest_qw2033.json`.",
        "2. Re-run in this order:",
        "   - `python3 QW_2030_FINAL_STAGE_C_GATE_COMBINED_BRANCH.py`",
        "   - `python3 QW_2031_V2_ETA_TRIAD_BLIND_EXTERNAL_VALIDATION.py`",
        "   - `python3 QW_2032_COMBINED_BRANCH_CONFIRMATORY_GATE.py`",
        "3. Confirm final verdict in `report_qw2032_combined_branch_confirmatory_gate.json`.",
        "",
        "## Expected Final Status",
        "- verdict: `COMBINED_BRANCH_CONFIRMATORY_GATE_PASS_STRONG`",
        "- readiness: `STAGE_C_PLUS_EXTERNAL_PRECONFIRMATORY_CLOSED`",
        "",
        "## External Data Sources (Not Frozen In Git)",
        "- Large raw archives are external by policy (no binary payload push in the bundle).",
        f"- Canonical source list: `{DATA_SOURCES_DOC}`",
        "",
        "### Required public sources",
    ]
    for src in external_sources:
        lines.extend(
            [
                f"- {src['name']}:",
                f"  - dataset: `{src['dataset']}`",
                f"  - url: `{src['url']}`",
                f"  - local example: `{src['local_example_path']}`",
                f"  - note: {src['note']}",
            ]
        )
    lines.extend(
        [
            "",
        "## Frozen File List",
        ]
    )
    for f in files:
        status = "OK" if f["exists"] else "MISSING"
        lines.append(f"- [{status}] `{f['path']}`")
    runbook.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return runbook


def main() -> None:
    BUNDLE_DIR.mkdir(parents=True, exist_ok=True)

    tracked_rel_paths = [
        "QW_2015_TRUE_EXTERNAL_BETA_CHANNEL_V2_READINESS_GATE.py",
        "QW_2017_V2_BETA_OBSERVABLE_BLIND_EXTERNAL_INTERVENTION.py",
        "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py",
        "QW_2027_GW_CONTROL_GAP_STRUCTURAL_TERM_SCAN.py",
        "QW_2029_CKM_BLOCKER_SHARED_FLAVOR_BASIS_SCAN.py",
        "QW_2030_FINAL_STAGE_C_GATE_COMBINED_BRANCH.py",
        "QW_2031_V2_ETA_TRIAD_BLIND_EXTERNAL_VALIDATION.py",
        "QW_2032_COMBINED_BRANCH_CONFIRMATORY_GATE.py",
        "report_qw2015_true_external_beta_channel_v2_readiness_gate.json",
        "report_qw2017_v2_beta_observable_blind_external_intervention.json",
        "report_qw2021_v2_eta_operator_beta_constraint_scan.json",
        "report_qw2027_gw_control_gap_structural_term_scan.json",
        "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json",
        "report_qw2030_final_stage_c_gate_combined_branch.json",
        "report_qw2031_v2_eta_triad_blind_external_validation.json",
        "report_qw2032_combined_branch_confirmatory_gate.json",
        "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv",
        "external_confirmatory_v2/beta_channel_true_external_v2/intervention_events.csv",
        "external_confirmatory_v2/beta_channel_true_external_v2/manifest_beta_channel.json",
        "external_confirmatory_v2/beta_channel_true_external_v2/protocol_freeze.json",
        "gw1831_window_features.csv",
    ]

    entries = [collect_file_entry(ROOT / rel) for rel in tracked_rel_paths]
    missing = [x["path"] for x in entries if not x["exists"]]

    runbook_path = write_runbook(BUNDLE_DIR, entries, EXTERNAL_REQUIRED_SOURCES)
    manifest_path = BUNDLE_DIR / "manifest_qw2033.json"

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "bundle_id": "QW2033_INDEPENDENT_CONFIRMATORY_FREEZE",
        "source_root": str(ROOT),
        "files": entries,
        "external_required_sources": EXTERNAL_REQUIRED_SOURCES,
        "external_sources_document": DATA_SOURCES_DOC,
        "missing_files": missing,
        "runbook": str(runbook_path.relative_to(ROOT)),
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    flags = {
        "all_required_files_present": bool(len(missing) == 0),
        "data_sources_doc_present": bool((ROOT / DATA_SOURCES_DOC).exists()),
        "bundle_manifest_written": bool(manifest_path.exists()),
        "bundle_runbook_written": bool(runbook_path.exists()),
    }
    pass_count = int(sum(1 for v in flags.values() if v))

    if pass_count == len(flags):
        verdict = "INDEPENDENT_CONFIRMATORY_FREEZE_BUNDLE_READY"
        readiness = "READY_FOR_TRUE_INDEPENDENT_MULTITEAM_RUN"
    else:
        verdict = "INDEPENDENT_CONFIRMATORY_FREEZE_BUNDLE_INCOMPLETE"
        readiness = "NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "bundle_dir": str(BUNDLE_DIR.relative_to(ROOT)),
        "manifest_path": str(manifest_path.relative_to(ROOT)),
        "runbook_path": str(runbook_path.relative_to(ROOT)),
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": len(flags),
        "missing_files": missing,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "HAND_OVER_BUNDLE_TO_INDEPENDENT_TEAM_FOR_BLIND_REPLICATION"
            if verdict == "INDEPENDENT_CONFIRMATORY_FREEZE_BUNDLE_READY"
            else "FILL_MISSING_FILES_AND_REBUILD_BUNDLE"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2033: INDEPENDENT CONFIRMATORY FREEZE BUNDLE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{len(flags)}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Bundle Artifacts",
            f"- bundle_dir: `{out['bundle_dir']}`",
            f"- manifest: `{out['manifest_path']}`",
            f"- runbook: `{out['runbook_path']}`",
            "",
            "## Missing Files",
        ]
    )
    if missing:
        for m in missing:
            lines.append(f"- {m}")
    else:
        lines.append("- none")
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

    print(f"[QW-2033] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2033] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2033] readiness={readiness} verdict={verdict} pass={pass_count}/{len(flags)}")


if __name__ == "__main__":
    main()
