#!/usr/bin/env python3
"""
QW-1900: Frozen external execution for QW-1899 protocol.

If external dataset package is not yet provided, emits a formal BLOCKED status.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
PROTO_PATH = ROOT / "external_confirmatory_v2" / "protocol_qw1899_external_detector.json"
MANIFEST_PATH = ROOT / "external_confirmatory_v2" / "manifest.json"
OUT_JSON = ROOT / "report_qw1900_frozen_external_run.json"
OUT_MD = ROOT / "RAPORT_QW1900_FROZEN_EXTERNAL_RUN.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    if not PROTO_PATH.exists():
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "FROZEN_EXTERNAL_RUN_BLOCKED",
            "reason": "protocol_missing",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1900: FROZEN EXTERNAL RUN\n\n"
            f"- Data UTC: {out['generated_utc']}\n"
            f"- Verdict: **{out['verdict']}**\n"
            "- Reason: protocol_qw1899_external_detector.json missing\n",
            encoding="utf-8",
        )
        print(f"[QW-1900] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1900] Saved MD:   {OUT_MD.name}")
        return

    proto_hash = sha256_file(PROTO_PATH)
    protocol = json.loads(PROTO_PATH.read_text(encoding="utf-8"))

    if not MANIFEST_PATH.exists():
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_sha256": proto_hash,
            "verdict": "FROZEN_EXTERNAL_RUN_BLOCKED",
            "readiness": "WAITING_EXTERNAL_DATASET_PACKAGE",
            "missing": ["external_confirmatory_v2/manifest.json", "pta_v2_pairs.csv", "gw_windows.csv"],
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

        lines = [
            "# RAPORT QW-1900: FROZEN EXTERNAL RUN",
            "",
            f"- Data UTC: {out['generated_utc']}",
            f"- Verdict: **{out['verdict']}**",
            f"- Readiness: **{out['readiness']}**",
            f"- Protocol SHA256: `{proto_hash}`",
            "",
            "## Missing Package",
            "- external_confirmatory_v2/manifest.json",
            "- external_confirmatory_v2/pta_v2_pairs.csv",
            "- external_confirmatory_v2/gw_windows.csv",
            "",
            "## Locked Rule",
            "- No external run before dataset freeze with manifest hashes.",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
        OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

        print(f"[QW-1900] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1900] Saved MD:   {OUT_MD.name}")
        return

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": proto_hash,
        "manifest_present": True,
        "verdict": "FROZEN_EXTERNAL_RUN_READY_TO_EXECUTE",
        "note": "Implement metric execution after data freeze validation.",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1900: FROZEN EXTERNAL RUN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- Protocol SHA256: `{proto_hash}`",
        "",
        "## Status",
        "- Manifest found. Proceed to frozen external metric execution (QW-1901).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1900] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1900] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
