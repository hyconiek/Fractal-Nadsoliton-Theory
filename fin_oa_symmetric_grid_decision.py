#!/usr/bin/env python3
"""Symmetric finite-grid composite decision with abstention."""
from __future__ import annotations

import math
from typing import Any

from fin_oa_mixture_eprocess import logsumexp


def score_symmetric(grid: dict[str, Any], events: list[dict[str, Any]], alpha: float=.01):
    components=grid["components"];pC=grid["p_C"]
    logs=[0.0 for _ in components];logweights=[math.log(float(c["weight"])) for c in components]
    upper=math.log(1/alpha);lower=math.log(alpha)
    for n,event in enumerate(events,1):
        ti=event.get("time_index");ret=event.get("return")
        if not isinstance(ti,int) or not 0<=ti<len(pC) or ret not in (0,1,False,True):
            return {"decision":"INVALID","n":n-1,"mixture_logLR":float('-inf'),"max_component_logLR":float('inf')}
        x=int(ret);pc=float(pC[ti])
        for j,c in enumerate(components):
            pq=float(c["p_Q"][ti]);logs[j]+=math.log(pq/pc) if x else math.log((1-pq)/(1-pc))
        mixture=logsumexp([lw+lr for lw,lr in zip(logweights,logs)])
        maximum=max(logs)
        if mixture>=upper:return {"decision":"Q_MIXTURE","n":n,"mixture_logLR":mixture,"max_component_logLR":maximum}
        if maximum<=lower:return {"decision":"C_ALL_GRID","n":n,"mixture_logLR":mixture,"max_component_logLR":maximum}
    return {"decision":"CONTINUE","n":len(events),"mixture_logLR":logsumexp([lw+lr for lw,lr in zip(logweights,logs)]),"max_component_logLR":max(logs)}

