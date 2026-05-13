#!/usr/bin/env python3
"""P1434: apply non-destructive pass-scope semantics retrofit to P1412-P1431 summaries."""

from __future__ import annotations

import json
from pathlib import Path

TARGET_PREFIXES = [f"p{n}" for n in range(1412, 1432)]


def main() -> None:
    root = Path(__file__).resolve().parent
    gen = root / "generated"
    updated = []

    for path in sorted(gen.glob("*.json")):
        name = path.name.lower()
        if not any(name.startswith(prefix) for prefix in TARGET_PREFIXES):
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        changed = False
        if "scope_of_pass" not in data:
            data["scope_of_pass"] = "local_contract"
            changed = True
        if "strict_core_qw2191_closed" not in data:
            data["strict_core_qw2191_closed"] = False
            changed = True
        if changed:
            path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
            updated.append(path.name)

    summary = {
        "packet": "P1434",
        "status": "PASS_SCOPE_RETROFIT_APPLIED",
        "route": "F-Nadsoliton=>L_SM+L_GR",
        "updated_files_count": len(updated),
        "updated_files": updated,
        "retrofit_values": {
            "scope_of_pass": "local_contract",
            "strict_core_qw2191_closed": False,
        },
        "legacy_import_used": False,
        "next_honest_step": "enforce these fields in all new checkpoints by default",
    }

    out = gen / "p1434_pass_scope_retrofit_apply_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[P1434] updated={len(updated)} wrote {out}")


if __name__ == "__main__":
    main()
