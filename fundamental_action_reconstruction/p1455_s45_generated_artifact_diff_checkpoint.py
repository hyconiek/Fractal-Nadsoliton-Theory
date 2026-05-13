#!/usr/bin/env python3
"""P1455 S4.5: compact local-only diff for generated JSON/CSV artifacts vs P1454."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
BASELINE = GEN / "p1454_generated_artifact_manifest.json"
CURRENT = GEN / "p1455_generated_artifact_manifest.json"
DIFF_JSON = GEN / "p1455_s45_generated_artifact_diff_summary.json"
DIFF_CSV = GEN / "p1455_s45_generated_artifact_diff_rows.csv"

EXCLUDED_PREFIXES = ("p1455_",)
EXCLUDED_EXACT = {
    "p1454_generated_artifact_manifest.json",
    "p1454_generated_artifact_manifest.csv",
}


def should_include(name: str, suffix: str) -> bool:
    if suffix not in {".json", ".csv"}:
        return False
    if name in EXCLUDED_EXACT:
        return False
    if any(name.startswith(prefix) for prefix in EXCLUDED_PREFIXES):
        return False
    return True


def collect_current_rows() -> list[dict]:
    rows: list[dict] = []
    for path in sorted(GEN.iterdir()):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        name = path.name
        if not should_include(name, suffix):
            continue
        rows.append({"name": name, "ext": suffix, "size_bytes": int(path.stat().st_size)})
    return rows


def normalize_rows(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for row in rows:
        name = row.get("name", row.get("file"))
        size = row.get("size_bytes", row.get("bytes"))
        ext = row.get("ext", Path(str(name)).suffix.lower() if name else None)
        if name is None or size is None:
            continue
        if not should_include(str(name), str(ext).lower()):
            continue
        out.append({"name": str(name), "ext": str(ext).lower(), "size_bytes": int(size)})
    return out


def build_digest(rows: list[dict]) -> str:
    canonical = "\n".join(f"{r['name']}|{r['size_bytes']}" for r in sorted(rows, key=lambda x: x["name"]))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def main() -> None:
    baseline_raw = json.loads(BASELINE.read_text(encoding="utf-8"))
    baseline_rows = normalize_rows(baseline_raw.get("artifacts", []))
    current_rows = collect_current_rows()

    b = {r["name"]: r for r in baseline_rows}
    c = {r["name"]: r for r in current_rows}

    added = sorted(set(c) - set(b))
    removed = sorted(set(b) - set(c))
    common = sorted(set(b) & set(c))

    changed: list[dict] = []
    unchanged = 0
    for name in common:
        before = int(b[name]["size_bytes"])
        after = int(c[name]["size_bytes"])
        if before != after:
            changed.append({"name": name, "size_before": before, "size_after": after, "delta": after - before})
        else:
            unchanged += 1

    status = "PASS_DIFF_STABLE" if not added and not removed and not changed else "PASS_DIFF_WITH_CHANGES"

    current_manifest = {
        "packet": "P1455",
        "scope": "generated json/csv only",
        "strict_only": True,
        "legacy_bridge_used": False,
        "excluded_exact": sorted(EXCLUDED_EXACT),
        "excluded_prefixes": sorted(EXCLUDED_PREFIXES),
        "total": len(current_rows),
        "digest_sha256": build_digest(current_rows),
        "sample_head": [r["name"] for r in sorted(current_rows, key=lambda x: x["name"])[:25]],
    }
    CURRENT.write_text(json.dumps(current_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    diff_summary = {
        "packet": "P1455",
        "status": status,
        "baseline_packet": "P1454",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "baseline_total": len(baseline_rows),
        "current_total": len(current_rows),
        "added_count": len(added),
        "removed_count": len(removed),
        "size_changed_count": len(changed),
        "unchanged_count": unchanged,
        "added": added,
        "removed": removed,
        "size_changed": changed,
    }
    DIFF_JSON.write_text(json.dumps(diff_summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    with DIFF_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["change_type", "name", "size_before", "size_after", "delta"])
        writer.writeheader()
        for name in added:
            writer.writerow({"change_type": "added", "name": name, "size_before": "", "size_after": c[name]["size_bytes"], "delta": ""})
        for name in removed:
            writer.writerow({"change_type": "removed", "name": name, "size_before": b[name]["size_bytes"], "size_after": "", "delta": ""})
        for row in changed:
            writer.writerow({"change_type": "size_changed", **row})

    print(f"[P1455] status={status} added={len(added)} removed={len(removed)} changed={len(changed)}")


if __name__ == "__main__":
    main()
