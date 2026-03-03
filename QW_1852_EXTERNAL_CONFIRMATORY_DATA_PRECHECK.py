#!/usr/bin/env python3
"""
QW-1852: External confirmatory data precheck package for JOINT V2.

Creates a frozen intake package (templates + README) and validates candidate
external dataset against protocol hashes, schema, and anti-leakage guards.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parent
PACKAGE_DIR = ROOT / "external_confirmatory_v2"
OUT_JSON = ROOT / "report_qw1852_external_confirmatory_data_precheck.json"
OUT_MD = ROOT / "RAPORT_QW1852_EXTERNAL_CONFIRMATORY_DATA_PRECHECK.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def write_csv_header(path: Path, header: List[str]) -> None:
    if path.exists():
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)


def ensure_package(pta_hash: str, gw_hash: str) -> Dict[str, str]:
    PACKAGE_DIR.mkdir(parents=True, exist_ok=True)

    readme_path = PACKAGE_DIR / "README_EXTERNAL_CONFIRMATORY_V2.md"
    manifest_template_path = PACKAGE_DIR / "manifest_template.json"
    pta_template_path = PACKAGE_DIR / "pta_v2_pairs_template.csv"
    gw_template_path = PACKAGE_DIR / "gw_windows_template.csv"

    if not readme_path.exists():
        readme_lines = [
            "# External Confirmatory V2 Intake Package",
            "",
            "This directory contains templates for a NEW external confirmatory dataset.",
            "Do not reuse design/training datasets.",
            "",
            "## Required Files",
            "- `manifest.json` copied from `manifest_template.json` and filled in.",
            "- `pta_v2_pairs.csv` with pair-level PTA features.",
            "- `gw_windows.csv` with GW window-level features.",
            "",
            "## Required Protocol Hashes",
            f"- PTA V2 protocol hash: `{pta_hash}`",
            f"- GW protocol hash: `{gw_hash}`",
            "",
            "## Minimal Rules",
            "- Data must be external to design dataset used in QW-1848/1850.",
            "- Keep file hashes stable after manifest freeze.",
            "- No post-hoc threshold changes.",
        ]
        readme_path.write_text("\n".join(readme_lines) + "\n", encoding="utf-8")

    if not manifest_template_path.exists():
        manifest_template = {
            "protocol": {
                "pta_v2_protocol_sha256": pta_hash,
                "gw_protocol_sha256": gw_hash,
            },
            "dataset": {
                "dataset_id": "EXTERNAL_DATASET_ID",
                "provider": "EXTERNAL_PROVIDER",
                "externality_statement": "Dataset is independent from QW-1848 design data.",
                "license": "TBD",
                "prepared_utc": "YYYY-MM-DDTHH:MM:SSZ",
            },
            "files": [
                {
                    "role": "pta_pairs",
                    "path": "pta_v2_pairs.csv",
                    "sha256": "FILL_AFTER_FILE_FREEZE",
                },
                {
                    "role": "gw_windows",
                    "path": "gw_windows.csv",
                    "sha256": "FILL_AFTER_FILE_FREEZE",
                },
            ],
        }
        manifest_template_path.write_text(json.dumps(manifest_template, ensure_ascii=False, indent=2), encoding="utf-8")

    write_csv_header(
        pta_template_path,
        [
            "pair_id",
            "theta_deg",
            "hxy",
            "f_mean",
            "f_std",
            "f_slope",
            "f_quad",
            "f_spread",
            "f_autoc1",
            "f_switch",
        ],
    )
    write_csv_header(
        gw_template_path,
        [
            "pair",
            "window_idx",
            "max_abs_corr",
            "mean_abs_corr",
            "corr_at_0ms",
            "corr_at_10ms",
        ],
    )

    return {
        "readme": str(readme_path.name),
        "manifest_template": str(manifest_template_path.name),
        "pta_template": str(pta_template_path.name),
        "gw_template": str(gw_template_path.name),
    }


def expected_design_hashes() -> Dict[str, str]:
    refs: Dict[str, str] = {}

    design_files = [
        ROOT / "gw1831_window_features.csv",
        ROOT / "NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz",
    ]
    for p in design_files:
        if p.exists() and p.is_file():
            refs[str(p.name)] = sha256_file(p)
    return refs


def validate_candidate(
    candidate_dir: Path,
    pta_hash: str,
    gw_hash: str,
    design_hashes: Dict[str, str],
) -> Dict:
    res = {
        "candidate_dir": str(candidate_dir),
        "exists": candidate_dir.exists(),
        "errors": [],
        "warnings": [],
        "checks": {},
    }

    if not candidate_dir.exists():
        res["errors"].append("candidate_dir_missing")
        return res

    manifest_path = candidate_dir / "manifest.json"
    if not manifest_path.exists():
        res["errors"].append("manifest_missing")
        return res

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    res["manifest"] = manifest

    m_pta = manifest.get("protocol", {}).get("pta_v2_protocol_sha256")
    m_gw = manifest.get("protocol", {}).get("gw_protocol_sha256")
    res["checks"]["protocol_hash_match"] = bool(m_pta == pta_hash and m_gw == gw_hash)
    if not res["checks"]["protocol_hash_match"]:
        res["errors"].append("protocol_hash_mismatch")

    ext_statement = str(manifest.get("dataset", {}).get("externality_statement", "")).strip()
    res["checks"]["externality_statement_present"] = len(ext_statement) >= 20
    if not res["checks"]["externality_statement_present"]:
        res["errors"].append("externality_statement_missing_or_too_short")

    files = manifest.get("files", [])
    role_map = {x.get("role"): x for x in files if isinstance(x, dict)}

    required_roles = {"pta_pairs", "gw_windows"}
    missing_roles = [r for r in required_roles if r not in role_map]
    if missing_roles:
        res["errors"].append(f"missing_roles:{','.join(sorted(missing_roles))}")
        return res

    seen_hashes = set()
    for role in sorted(required_roles):
        entry = role_map[role]
        rel_path = entry.get("path")
        expected_hash = str(entry.get("sha256", ""))
        if not rel_path:
            res["errors"].append(f"{role}:path_missing")
            continue

        fpath = (candidate_dir / rel_path).resolve()
        try:
            fpath.relative_to(candidate_dir.resolve())
        except Exception:
            res["errors"].append(f"{role}:path_outside_candidate_dir")
            continue

        if not fpath.exists():
            res["errors"].append(f"{role}:file_missing")
            continue

        real_hash = sha256_file(fpath)
        seen_hashes.add(real_hash)
        hash_ok = real_hash == expected_hash
        res["checks"][f"{role}_hash_match"] = bool(hash_ok)
        if not hash_ok:
            res["errors"].append(f"{role}:sha256_mismatch")

        if real_hash in design_hashes.values():
            res["errors"].append(f"{role}:matches_design_dataset_hash")

        df = pd.read_csv(fpath)
        if role == "pta_pairs":
            req = [
                "theta_deg",
                "hxy",
                "f_mean",
                "f_std",
                "f_slope",
                "f_quad",
                "f_spread",
                "f_autoc1",
                "f_switch",
            ]
            missing = [c for c in req if c not in df.columns]
            if missing:
                res["errors"].append("pta_pairs:missing_columns:" + ",".join(missing))
            else:
                if len(df) < 80:
                    res["errors"].append("pta_pairs:too_few_rows_min80")
                if df[req].isnull().any().any():
                    res["errors"].append("pta_pairs:contains_nan")
                res["checks"]["pta_pairs_row_count"] = int(len(df))

        elif role == "gw_windows":
            req = [
                "pair",
                "window_idx",
                "max_abs_corr",
                "mean_abs_corr",
                "corr_at_0ms",
                "corr_at_10ms",
            ]
            missing = [c for c in req if c not in df.columns]
            if missing:
                res["errors"].append("gw_windows:missing_columns:" + ",".join(missing))
            else:
                if len(df) < 150:
                    res["errors"].append("gw_windows:too_few_rows_min150")
                if df[req].isnull().any().any():
                    res["errors"].append("gw_windows:contains_nan")
                pairs = set(df["pair"].astype(str).unique().tolist())
                required_pairs = {"H1-L1", "H1-V1", "L1-V1"}
                if not required_pairs.issubset(pairs):
                    res["errors"].append("gw_windows:missing_required_detector_pairs")
                res["checks"]["gw_windows_row_count"] = int(len(df))

    if len(seen_hashes) != len(required_roles):
        res["warnings"].append("file_hash_collision_or_missing_hashes")

    res["ok"] = len(res["errors"]) == 0
    return res


def main() -> None:
    parser = argparse.ArgumentParser(description="QW-1852 external confirmatory data precheck")
    parser.add_argument(
        "--candidate-dir",
        type=str,
        default=str(PACKAGE_DIR / "confirmatory_dataset"),
        help="Directory containing manifest.json and data files for external confirmatory run.",
    )
    args = parser.parse_args()

    d1839 = load_json("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1850 = load_json("report_qw1850_pta_v2_prereg_protocol.json")

    gw_hash = str(d1839.get("protocol_sha256"))
    pta_hash = str(d1850.get("protocol_sha256"))

    package_files = ensure_package(pta_hash=pta_hash, gw_hash=gw_hash)
    design_hashes = expected_design_hashes()

    candidate_dir = Path(args.candidate_dir)
    validation = validate_candidate(candidate_dir, pta_hash=pta_hash, gw_hash=gw_hash, design_hashes=design_hashes)

    if validation.get("ok"):
        readiness = "EXTERNAL_DATASET_READY_FOR_QW1853"
        hard_gate = "PASS"
    elif validation.get("exists") and ("manifest_missing" not in validation.get("errors", [])):
        readiness = "EXTERNAL_DATASET_PRESENT_BUT_INVALID"
        hard_gate = "FAIL"
    else:
        readiness = "EXTERNAL_DATASET_PENDING_COLLECTION"
        hard_gate = "PARTIAL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_hashes": {
            "pta_v2_protocol_sha256": pta_hash,
            "gw_protocol_sha256": gw_hash,
        },
        "package_dir": str(PACKAGE_DIR),
        "package_files": package_files,
        "design_reference_hashes": design_hashes,
        "validation": validation,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "verdict": "EXTERNAL_CONFIRMATORY_DATA_PRECHECK_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1852: EXTERNAL CONFIRMATORY DATA PRECHECK",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Protocol Hashes",
        f"- PTA V2: `{pta_hash}`",
        f"- GW: `{gw_hash}`",
        "",
        "## Package",
        f"- Directory: `{PACKAGE_DIR}`",
        f"- Files: {', '.join(package_files.values())}",
        "",
        "## Candidate Validation",
        f"- Candidate dir: `{validation.get('candidate_dir')}`",
        f"- Exists: {validation.get('exists')}",
        f"- OK: {validation.get('ok', False)}",
        f"- Errors: {validation.get('errors', [])}",
        f"- Warnings: {validation.get('warnings', [])}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1852] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1852] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
