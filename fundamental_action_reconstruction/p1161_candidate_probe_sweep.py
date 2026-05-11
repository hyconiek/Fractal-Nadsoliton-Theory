#!/usr/bin/env python3
"""P1161: physically-motivated candidate sweep under strict-side hints.

Computes P1146-like observables per candidate (without forcing stage failure)
and ranks by lower sign-change burden then lower negative_count.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
ALPHA = 4.0 * math.log(2.0)


def eval_candidate(payload: dict) -> dict:
    omega = float(payload.get("omega_hint", 0.18575))
    phi = float(payload.get("phi_hint", 0.16250))
    beta = float(payload.get("beta_hint", 1.0))
    eta = float(payload.get("eta_hint", 1.8))

    vals=[]
    for i in range(241):
        d=i*0.1
        v=math.exp(-ALPHA*d)*math.cos(omega*d+phi)/(1.0+beta*(d**eta))
        vals.append(v)

    pos=sum(v>0 for v in vals); neg=sum(v<0 for v in vals)
    sign_changes=0; prev=0
    for v in vals:
        s=1 if v>0 else (-1 if v<0 else 0)
        if s!=0 and prev!=0 and s!=prev: sign_changes+=1
        if s!=0: prev=s
    return {
        "omega": omega, "phi": phi, "beta": beta, "eta": eta,
        "negative_count": neg,
        "positive_count": pos,
        "sign_change_count": sign_changes,
        "qw_2191_proxy": "BLOCKED" if (neg>0 or sign_changes>0) else "CANDIDATE_PASS_ONLY",
    }


def main() -> None:
    reg = json.loads((GEN / "p1161_candidate_registry_physical.json").read_text())
    rows=[]
    for c in reg["candidates"]:
        p=json.loads(Path(c).read_text())
        r=eval_candidate(p)
        r["candidate"]=c
        r["premise_id"]=p.get("premise_id")
        rows.append(r)
    rows.sort(key=lambda x:(x["sign_change_count"], x["negative_count"]))
    out={
      "packet":"P1161","as_of":"2026-05-10","rows":rows,
      "top_recommendation": rows[0]["candidate"] if rows else None,
      "note":"Candidate investigation only; no closure/QW-2191 discharge claim."
    }
    outp=GEN/"p1161_candidate_probe_sweep_summary.json"
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(f"[P1161] evaluated {len(rows)} candidates wrote {outp}")

if __name__=='__main__':
    main()
