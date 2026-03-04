#!/usr/bin/env python3
"""
QW-2054: Protocol lock integrity gate (tamper detection).

Checks:
- protocol_lock_qw2053.json canonical hash consistency,
- artifact hash list consistency against current files,
- runbook presence for external handoff.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
LOCK_JSON = ROOT / "external_confirmatory_v2" / "independent_multiteam_lock_qw2053" / "protocol_lock_qw2053.json"
LOCK_RUNBOOK = ROOT / "external_confirmatory_v2" / "independent_multiteam_lock_qw2053" / "RUNBOOK_QW2053.md"

OUT_JSON = ROOT / "report_qw2054_protocol_lock_integrity_gate.json"
OUT_MD = ROOT / "RAPORT_QW2054_PROTOCOL_LOCK_INTEGRITY_GATE.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha256(data: Dict) -> str:
    b = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(b).hexdigest()


def main() -> None:
    if not LOCK_JSON.exists():
        raise RuntimeError(f"Missing lock file: {LOCK_JSON}")

    lock = json.loads(LOCK_JSON.read_text(encoding="utf-8"))
    stored_lock_hash = str(lock.get("lock_sha256", ""))

    body = {k: v for k, v in lock.items() if k != "lock_sha256"}
    recomputed_lock_hash = canonical_sha256(body)

    artifact_rows: List[Dict] = lock.get("artifact_hashes", [])
    missing_files: List[str] = []
    hash_mismatches: List[Dict[str, str]] = []

    for row in artifact_rows:
        rel = str(row.get("path", ""))
        want = str(row.get("sha256", ""))
        p = ROOT / rel
        if not p.exists():
            missing_files.append(rel)
            continue
        got = sha256_file(p)
        if got != want:
            hash_mismatches.append({"path": rel, "expected": want, "got": got})

    flags = {
        "lock_file_exists": LOCK_JSON.exists(),
        "runbook_exists": LOCK_RUNBOOK.exists(),
        "stored_lock_hash_matches_recomputed": stored_lock_hash == recomputed_lock_hash,
        "all_artifact_files_present": len(missing_files) == 0,
        "all_artifact_hashes_match": len(hash_mismatches) == 0,
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if pass_count == total_flags:
        verdict = "PROTOCOL_LOCK_INTEGRITY_PASS"
        readiness = "LOCK_IS_TAMPER_EVIDENT_AND_READY"
    elif pass_count >= total_flags - 1:
        verdict = "PROTOCOL_LOCK_INTEGRITY_PARTIAL"
        readiness = "PARTIAL"
    else:
        verdict = "PROTOCOL_LOCK_INTEGRITY_FAIL"
        readiness = "NOT_READY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "lock_file": str(LOCK_JSON.relative_to(ROOT)),
        "stored_lock_sha256": stored_lock_hash,
        "recomputed_lock_sha256": recomputed_lock_hash,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "missing_files": missing_files,
        "hash_mismatches": hash_mismatches,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": "EXECUTE_TRULY_INDEPENDENT_MULTITEAM_CONFIRMATORY_RUN_UNDER_LOCK_QW2053",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2054: PROTOCOL LOCK INTEGRITY GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Lock Hash",
        f"- stored: `{stored_lock_hash}`",
        f"- recomputed: `{recomputed_lock_hash}`",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")

    if missing_files:
        lines.append("")
        lines.append("## Missing Files")
        for x in missing_files:
            lines.append(f"- {x}")

    if hash_mismatches:
        lines.append("")
        lines.append("## Hash Mismatches")
        for x in hash_mismatches:
            lines.append(f"- `{x['path']}`")

    lines.extend(
        [
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2054] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2054] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2054] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
