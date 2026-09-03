#!/usr/bin/env python3
"""Dependency-light fail-closed validator/scorer for FIN OA protocol 10.98."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


REQUIRED_PROTOCOL = {
    "schema_version", "matrix_fingerprint_sha256", "primary", "robust",
    "raw_event_fields", "invalid_record_rules", "role_requirements",
}


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def validate_protocol(protocol: dict[str, Any]) -> list[str]:
    errors=[]
    missing=sorted(REQUIRED_PROTOCOL-set(protocol))
    if missing: errors.append(f"missing protocol fields: {missing}")
    primary=protocol.get("primary",{})
    for key in ["time","shots","return_threshold_Q"]:
        if key not in primary: errors.append(f"missing primary.{key}")
    robust=protocol.get("robust",{})
    times=robust.get("times",[]);shots=robust.get("shots_each",[])
    if len(times)!=4: errors.append("robust.times must have length 4")
    if len(shots)!=4 or any(int(n)<=0 for n in shots): errors.append("robust.shots_each must contain four positive counts")
    if protocol.get("role_requirements",{}).get("provider_equals_analyst_allowed",True): errors.append("provider/analyst separation must be required")
    return errors


def validate_aggregate_record(protocol: dict[str, Any], record: dict[str, Any]) -> list[str]:
    errors=[]
    for key in ["model_blind_id","times","attempts","clicks","return_counts","config_id","run_id"]:
        if key not in record: errors.append(f"missing record field: {key}")
    times=record.get("times",[]);attempts=record.get("attempts",[]);clicks=record.get("clicks",[]);returns=record.get("return_counts",[])
    if not (len(times)==len(attempts)==len(clicks)==len(returns)): errors.append("record vector lengths differ")
    for i,(a,c,r) in enumerate(zip(attempts,clicks,returns)):
        if not all(isinstance(v,int) for v in (a,c,r)): errors.append(f"noninteger count at {i}");continue
        if a<0 or c<0 or r<0 or r>c or c>a: errors.append(f"invalid count ordering at {i}")
    return errors


def score_primary(protocol: dict[str, Any], record: dict[str, Any]) -> str:
    pe=validate_protocol(protocol);re=validate_aggregate_record(protocol,record)
    if pe or re: return "INVALID"
    target=float(protocol["primary"]["time"])
    try: idx=[float(t) for t in record["times"]].index(target)
    except ValueError: return "INVALID"
    if record["clicks"][idx] != int(protocol["primary"]["shots"]): return "INVALID"
    return "Q" if record["return_counts"][idx] >= int(protocol["primary"]["return_threshold_Q"]) else "C"

