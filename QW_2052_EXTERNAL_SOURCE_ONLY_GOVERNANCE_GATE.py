#!/usr/bin/env python3
"""
QW-2052: External-source-only governance gate.

Goal:
- enforce "no large binary payload in git freeze bundles",
- require explicit source-document references for external public datasets,
- keep reproducibility based on manifests + scripts + source provenance.
"""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2052_external_source_only_governance_gate.json"
OUT_MD = ROOT / "RAPORT_QW2052_EXTERNAL_SOURCE_ONLY_GOVERNANCE_GATE.md"
DATA_DOC = ROOT / "DATA_SOURCES_EXTERNAL_DOWNLOADS.md"

RUNBOOKS = [
    ROOT / "external_confirmatory_v2" / "independent_bundle_qw2033" / "RUNBOOK_QW2033.md",
    ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge" / "RUNBOOK_QW2050.md",
]

RELEASE_DOCS = [
    ROOT / "RELEASE_4_9.md",
    ROOT / "RELEASE_4_9_TEXTBOOK_EN_PL.md",
]

MANIFEST_QW2033 = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2033" / "manifest_qw2033.json"
MANIFEST_QW2050 = ROOT / "external_confirmatory_v2" / "independent_bundle_qw2050_spectral_micro_bridge" / "manifest_qw2050.json"

SIZE_THRESHOLD_BYTES = 10 * 1024 * 1024


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception:
        return ""


def tracked_large_files() -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    raw = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True)
    for rel in raw.splitlines():
        p = ROOT / rel
        if not p.exists() or not p.is_file():
            continue
        sz = int(p.stat().st_size)
        if sz >= SIZE_THRESHOLD_BYTES:
            out.append({"path": rel, "bytes": sz})
    out.sort(key=lambda x: int(x["bytes"]), reverse=True)
    return out


def large_files_in_bundle_dirs() -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    parent = ROOT / "external_confirmatory_v2"
    for p in sorted(parent.rglob("independent_bundle_qw*/*")):
        if not p.is_file():
            continue
        sz = int(p.stat().st_size)
        if sz >= SIZE_THRESHOLD_BYTES:
            out.append({"path": str(p.relative_to(ROOT)), "bytes": sz})
    return out


def contains_reference(path: Path, token: str) -> bool:
    return token in read_text(path)


def load_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def any_large_archive_in_manifest(file_entries: List[Dict]) -> bool:
    for row in file_entries:
        rel = str(row.get("path", ""))
        if rel.endswith(".tar.gz") or rel.endswith(".h5") or rel.endswith(".hdf5"):
            return True
    return False


def main() -> None:
    large_tracked = tracked_large_files()
    large_in_bundles = large_files_in_bundle_dirs()

    m2033 = load_json(MANIFEST_QW2033)
    m2050 = load_json(MANIFEST_QW2050)

    m2033_files = m2033.get("files", [])
    m2050_files = m2050.get("manifest", [])

    flags = {
        "data_sources_doc_exists": DATA_DOC.exists(),
        "runbooks_reference_data_sources_doc": all(
            rb.exists() and contains_reference(rb, DATA_DOC.name) for rb in RUNBOOKS
        ),
        "release_docs_reference_data_sources_doc": all(
            rd.exists() and contains_reference(rd, DATA_DOC.name) for rd in RELEASE_DOCS
        ),
        "manifest_qw2033_declares_external_sources": bool(m2033.get("external_required_sources")),
        "manifest_qw2050_declares_data_sources_doc": str(m2050.get("external_sources_document", "")) == DATA_DOC.name,
        "manifest_qw2033_has_no_embedded_large_archives": not any_large_archive_in_manifest(m2033_files),
        "manifest_qw2050_has_no_embedded_large_archives": not any_large_archive_in_manifest(m2050_files),
        "independent_bundle_dirs_have_no_large_payloads": len(large_in_bundles) == 0,
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "EXTERNAL_SOURCE_ONLY_GOVERNANCE_PASS"
        readiness = "SOURCE_ONLY_CONFIRMATORY_GOVERNANCE_READY"
    elif pass_count >= total_flags - 2:
        verdict = "EXTERNAL_SOURCE_ONLY_GOVERNANCE_PARTIAL"
        readiness = "PARTIAL"
    else:
        verdict = "EXTERNAL_SOURCE_ONLY_GOVERNANCE_FAIL"
        readiness = "NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "size_threshold_bytes": SIZE_THRESHOLD_BYTES,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "tracked_large_files": large_tracked,
        "large_files_in_independent_bundle_dirs": large_in_bundles,
        "tracked_large_files_count": len(large_tracked),
        "verdict": verdict,
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2052: EXTERNAL SOURCE-ONLY GOVERNANCE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- size_threshold_bytes: {SIZE_THRESHOLD_BYTES}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    lines.extend(
        [
            "",
            f"## Tracked Large Files (>= {SIZE_THRESHOLD_BYTES} B)",
            f"- count: {len(large_tracked)}",
        ]
    )
    if large_tracked:
        for row in large_tracked[:30]:
            lines.append(f"- `{row['path']}` ({row['bytes']} B)")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Independent Bundle Large Payload Check",
            f"- large files in independent bundle dirs: {len(large_in_bundles)}",
        ]
    )
    if large_in_bundles:
        for row in large_in_bundles:
            lines.append(f"- `{row['path']}` ({row['bytes']} B)")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Required Next Step",
            "- Keep freeze bundles source-only and provide this data-source document with each handoff.",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2052] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2052] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2052] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
