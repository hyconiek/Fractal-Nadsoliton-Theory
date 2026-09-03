#!/usr/bin/env python3
"""Fail-closed sequential validator/scorer for FIN OA protocol 10.99."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def validate_protocol(protocol: dict[str, Any]) -> list[str]:
    errors=[]
    for key in ["schema_version","decimal_spectrum","simple_full_record","sequential","composite","calibration","role_requirements"]:
        if key not in protocol:errors.append(f"missing {key}")
    seq=protocol.get("sequential",{})
    llr=seq.get("distance_class_log_likelihood_ratios",[])
    if len(llr)!=7:errors.append("sequential LLR vector must have length 7")
    alpha=seq.get("alpha")
    if not isinstance(alpha,(int,float)) or not (0<alpha<1):errors.append("alpha must lie in (0,1)")
    roles=protocol.get("role_requirements",{})
    if roles.get("provider_equals_analyst_allowed",True):errors.append("independent roles must be required")
    return errors


def validate_sequence(protocol: dict[str, Any], outcomes: list[Any]) -> list[str]:
    errors=validate_protocol(protocol)
    for i,x in enumerate(outcomes):
        if not isinstance(x,int) or not 0<=x<7:errors.append(f"invalid distance class at {i}")
    return errors


def score_sequence(protocol: dict[str, Any], outcomes: list[Any]) -> dict[str, Any]:
    errors=validate_sequence(protocol,outcomes)
    if errors:return {"decision":"INVALID","errors":errors,"n":0,"logLR":0.0}
    llr=protocol["sequential"]["distance_class_log_likelihood_ratios"]
    alpha=float(protocol["sequential"]["alpha"])
    upper=math.log(1/alpha);lower=math.log(alpha);total=0.0
    for n,x in enumerate(outcomes,1):
        total+=float(llr[x])
        if total>=upper:return {"decision":"Q","errors":[],"n":n,"logLR":total}
        if total<=lower:return {"decision":"C","errors":[],"n":n,"logLR":total}
    return {"decision":"CONTINUE","errors":[],"n":len(outcomes),"logLR":total}

