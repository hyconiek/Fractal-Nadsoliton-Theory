#!/usr/bin/env python3
"""P1435: enforce PASS-scope semantics gate for strict-route generated summaries."""

from __future__ import annotations

import json
from pathlib import Path

START = 1412


def extract_packet_number(filename: str) -> int | None:
    if not filename.startswith("p"):
        return None
    digits = []
    for ch in filename[1:]:
        if ch.isdigit():
            digits.append(ch)
        else:
            break
    if not digits:
        return None
    return int("".join(digits))


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    audited, violations = 0, []

    for path in sorted(gen.glob("p*.json")):
        if path.name == "p1435_pass_scope_enforcement_gate_summary.json":
            continue
        pkt = extract_packet_number(path.name)
        if pkt is None or pkt < START:
            continue
        audited += 1
        data = json.loads(path.read_text(encoding="utf-8"))
        miss = [k for k in ("scope_of_pass", "strict_core_qw2191_closed") if k not in data]
        if miss:
            violations.append({"file": path.name, "missing_fields": miss})

    status = "PASS_ENFORCEMENT_GATE" if not violations else "FAIL_ENFORCEMENT_GATE"
    summary = {
        "packet": "P1435",
        "status": status,
        "audited_file_count": audited,
        "violations_count": len(violations),
        "violations": violations,
        "required_fields": ["scope_of_pass", "strict_core_qw2191_closed"],
        "legacy_import_used": False,
        "route": "F-Nadsoliton=>L_SM+L_GR",
        "next_honest_step": "if PASS: keep gate as CI prerequisite; if FAIL: patch missing files before accepting new checkpoints",
        "scope_of_pass": "local_contract",
        "strict_core_qw2191_closed": False,
    }

    out = gen / "p1435_pass_scope_enforcement_gate_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[P1435] status={status} audited={audited} violations={len(violations)}")


if __name__ == "__main__":
    main()
