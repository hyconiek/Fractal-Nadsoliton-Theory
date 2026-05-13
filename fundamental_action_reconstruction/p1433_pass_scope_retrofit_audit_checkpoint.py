#!/usr/bin/env python3
"""P1433: audit retrofit readiness for PASS scope semantics in P1412-P1431 summaries."""

from __future__ import annotations

import json
from pathlib import Path

TARGET_PREFIXES = [f"p{n}" for n in range(1412, 1432)]


def main() -> None:
    root = Path(__file__).resolve().parent
    gen = root / "generated"
    missing = []
    audited = 0

    for path in sorted(gen.glob("*.json")):
        name = path.name.lower()
        if not any(name.startswith(prefix) for prefix in TARGET_PREFIXES):
            continue
        audited += 1
        data = json.loads(path.read_text(encoding="utf-8"))
        miss = []
        if "scope_of_pass" not in data:
            miss.append("scope_of_pass")
        if "strict_core_qw2191_closed" not in data:
            miss.append("strict_core_qw2191_closed")
        if miss:
            missing.append({"file": path.name, "missing_fields": miss})

    summary = {
        "packet": "P1433",
        "status": "RETROFIT_AUDIT_COMPLETE",
        "route": "F-Nadsoliton=>L_SM+L_GR",
        "audited_file_count": audited,
        "files_missing_semantics_fields": missing,
        "legacy_import_used": False,
        "next_honest_step": "apply non-destructive retrofit by adding scope_of_pass=local_contract and strict_core_qw2191_closed=false to all listed files",
    }

    out = gen / "p1433_pass_scope_retrofit_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[P1433] audited={audited} missing={len(missing)}")
    print(f"[P1433] wrote {out}")


if __name__ == "__main__":
    main()
