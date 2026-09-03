#!/usr/bin/env python3
"""Finite-mixture e-process for classical C versus predeclared dephased-Q grid."""
from __future__ import annotations

import math
from typing import Any


def logsumexp(values: list[float]) -> float:
    m=max(values)
    return m+math.log(sum(math.exp(v-m) for v in values))


def score_mixture(grid: dict[str, Any], events: list[dict[str, Any]], alpha: float=.01) -> dict[str, Any]:
    components=grid["components"];times=grid["times"];pC=grid["p_C"]
    logs=[math.log(float(c["weight"])) for c in components]
    threshold=math.log(1/alpha)
    for n,event in enumerate(events,1):
        ti=event.get("time_index");ret=event.get("return")
        if not isinstance(ti,int) or not 0<=ti<len(times) or ret not in (0,1,False,True):
            return {"decision":"INVALID","n":n-1,"log_e":float('-inf')}
        x=int(ret);pc=float(pC[ti])
        for j,component in enumerate(components):
            pq=float(component["p_Q"][ti])
            logs[j]+=math.log(pq/pc) if x else math.log((1-pq)/(1-pc))
        loge=logsumexp(logs)
        if loge>=threshold:return {"decision":"REJECT_C_FOR_MIXTURE_Q","n":n,"log_e":loge}
    return {"decision":"CONTINUE","n":len(events),"log_e":logsumexp(logs)}

