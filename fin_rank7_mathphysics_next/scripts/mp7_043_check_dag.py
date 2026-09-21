#!/usr/bin/env python3
import json
from pathlib import Path
R=Path(__file__).resolve().parents[1]
d=json.loads((R/'audit/MP7-043_dependency_graph.json').read_text())
nodes={n['id'] for n in d['nodes']}
adj={n:[] for n in nodes}; indeg={n:0 for n in nodes}
for a,b in d['edges']:
    if a not in nodes or b not in nodes: raise SystemExit(f'unknown node {a}->{b}')
    adj[a].append(b); indeg[b]+=1
q=sorted([n for n,v in indeg.items() if v==0]); order=[]
while q:
    n=q.pop(0); order.append(n)
    for m in adj[n]:
        indeg[m]-=1
        if indeg[m]==0: q.append(m); q.sort()
if len(order)!=len(nodes):
    cyc=sorted(n for n,v in indeg.items() if v>0)
    raise SystemExit('cycle detected: '+','.join(cyc))
print(json.dumps({'acyclic':True,'node_count':len(nodes),'edge_count':len(d['edges']),'topological_order':order},indent=2))
