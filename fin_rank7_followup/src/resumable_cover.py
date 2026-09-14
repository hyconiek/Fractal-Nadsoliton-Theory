from __future__ import annotations
import json, os, time
from pathlib import Path

def _leaf_value(a,b):
    # deterministic diagnostic polynomial; not a scientific claim
    m=(a+b)/2
    return (m-0.37)**2 + 0.01

def run(checkpoint:Path, max_leaves=128, stop_after=None):
    if checkpoint.exists():
        st=json.loads(checkpoint.read_text())
    else:
        st={'queue':[[0.0,1.0]],'leaves':[],'processed':0,'complete':False}
    while st['queue'] and len(st['leaves'])<max_leaves:
        a,b=st['queue'].pop(0); m=(a+b)/2
        # split until width <= 1/64, then accept diagnostic leaf
        if b-a>1/64:
            st['queue'].append([a,m]); st['queue'].append([m,b])
        else:
            st['leaves'].append({'a':a,'b':b,'value':_leaf_value(a,b)})
        st['processed']+=1
        tmp=checkpoint.with_suffix('.tmp')
        tmp.write_text(json.dumps(st,sort_keys=True))
        os.replace(tmp,checkpoint)
        if stop_after is not None and st['processed']>=stop_after:
            return st
    st['complete']=not st['queue']
    checkpoint.write_text(json.dumps(st,sort_keys=True))
    return st
