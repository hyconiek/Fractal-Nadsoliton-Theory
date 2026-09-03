#!/usr/bin/env python3
"""Fail-closed schema validator for the frozen ST276 child-resolved record.

This validates record structure and custody.  It never converts synthetic data
into empirical evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path


HEX64 = re.compile(r"^[0-9a-f]{64}$")


def validate(record: dict, require_complete: bool = True, empirical: bool = False) -> list[str]:
    errors: list[str] = []
    required = ["schema_version", "provider", "registrar", "analyst", "holdout_frozen",
                "calibrated_time", "run_id", "calibration_hash", "protocol_hash", "events"]
    for key in required:
        if key not in record:
            errors.append(f"missing field: {key}")
    if errors:
        return errors
    if record["schema_version"] != "FIN-ST276-1":
        errors.append("unsupported schema_version")
    roles = [record["provider"], record["registrar"], record["analyst"]]
    if len(set(roles)) != 3:
        errors.append("provider, registrar, and analyst must be distinct")
    if record["holdout_frozen"] is not True:
        errors.append("holdout is not frozen")
    if not isinstance(record["calibrated_time"], (int, float)) or record["calibrated_time"] <= 0:
        errors.append("calibrated_time must be positive")
    for key in ("calibration_hash", "protocol_hash"):
        if not isinstance(record[key], str) or not HEX64.fullmatch(record[key]):
            errors.append(f"{key} must be a lowercase SHA-256 digest")
    if empirical and record.get("synthetic", False):
        errors.append("synthetic record cannot pass the empirical-evidence gate")
    settings = set()
    for i, event in enumerate(record["events"]):
        keys = {"timestamp", "preparation_x", "preparation_child", "effect_y", "effect_child", "count"}
        if not keys.issubset(event):
            errors.append(f"event {i} has incomplete fields")
            continue
        x, a, y, b = (event["preparation_x"], event["preparation_child"],
                      event["effect_y"], event["effect_child"])
        if x not in range(12) or y not in range(12) or a not in (0, 1) or b not in (0, 1):
            errors.append(f"event {i} has an out-of-range setting")
        if not isinstance(event["count"], int) or event["count"] < 0:
            errors.append(f"event {i} has an invalid count")
        settings.add((x, a, y, b))
    if require_complete and len(settings) != 24 * 24:
        errors.append(f"incomplete child-resolved setting coverage: {len(settings)}/576")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("record", type=Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--empirical", action="store_true")
    args = parser.parse_args()
    record = json.loads(args.record.read_text(encoding="utf-8"))
    errors = validate(record, not args.allow_incomplete, args.empirical)
    payload = {"record": args.record.name, "sha256": hashlib.sha256(args.record.read_bytes()).hexdigest(),
               "valid": not errors, "errors": errors}
    print(json.dumps(payload, indent=2))
    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
