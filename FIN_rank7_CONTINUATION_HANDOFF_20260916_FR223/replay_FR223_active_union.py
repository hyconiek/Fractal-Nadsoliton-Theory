#!/usr/bin/env python3
"""Replay the active local union encoded in the FR223 residual-search artifact.

This is a portability helper, not a replacement for the original per-FR tests.
The FR223 JSON stores mask endpoints as JSON numbers, so this script replays the
exact decimal values serialized there (Fraction(str(value))). A PASS therefore
certifies the stored FR223 mask geometry itself.
"""
from __future__ import annotations
import json
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'src'
sys.path.insert(0,str(SRC))
import frontier_local_boxes as flb
import frontier_shifted_boxes as fsb

J=ROOT/'results'/'FR223_post_FR222_residual_search.json'
d=json.loads(J.read_text())

def q(x): return F(str(x))

rows=[]
for m in d['centered_masks']:
    out=flb.raw_box(q(m['rx']),q(m['ru']),q(m['rv']),q(m['e']))
    ok=out.get('status')=='INTERVAL_CERTIFIED'
    rows.append({'name':m['name'],'kind':'centered','ok':ok,'status':out.get('status'),'reason':out.get('reason')})

for m in d['shifted_masks']:
    box=((q(m['x'][0]),q(m['x'][1])),(q(m['u'][0]),q(m['u'][1])),(q(m['v'][0]),q(m['v'][1])),(F(0),q(m['e'])))
    out=fsb.raw_shifted_box(box)
    ok=out.get('status')=='INTERVAL_CERTIFIED'
    rows.append({'name':m['name'],'kind':'shifted','ok':ok,'status':out.get('status'),'reason':out.get('reason')})

bad=[r for r in rows if not r['ok']]
summary={'source':str(J.name),'centered_count':len(d['centered_masks']),'shifted_count':len(d['shifted_masks']),'total':len(rows),'pass':len(rows)-len(bad),'fail':len(bad),'failed':bad,'rows':rows}
outp=ROOT/'results'/'FR223_ACTIVE_UNION_REPLAY.json'
outp.write_text(json.dumps(summary,indent=2,sort_keys=True))
print(json.dumps({k:summary[k] for k in ['centered_count','shifted_count','total','pass','fail']},indent=2))
if bad:
    print('FAILED:', ', '.join(r['name'] for r in bad))
    raise SystemExit(1)
print('FR223_ACTIVE_UNION_REPLAY_PASS')
