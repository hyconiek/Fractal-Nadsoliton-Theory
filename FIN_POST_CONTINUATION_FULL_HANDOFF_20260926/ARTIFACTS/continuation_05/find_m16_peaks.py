#!/usr/bin/env python3
import subprocess,re,math
from pathlib import Path
BIN='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_05/fin12_m16_jet4'
ORBIT='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_03/d12_four_subset_orbits.tsv'
OUT='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_05/m16_peak_all29.tsv'
pat=re.compile(r'beta ([\d.eE+-]+) mask (\d+).*? C ([\d.eE+-]+) dC ([\d.eE+-]+).*? skew ([\d.eE+-]+) binder ([\d.eE+-]+)')
cache={}
def evalp(mask,b):
    k=(mask,round(b,12))
    if k in cache:return cache[k]
    s=subprocess.check_output([BIN,f'{b:.12g}',str(mask)],text=True)
    m=pat.search(s)
    if not m: raise RuntimeError(s)
    r=dict(beta=float(m.group(1)),C=float(m.group(3)),dC=float(m.group(4)),skew=float(m.group(5)),binder=float(m.group(6)))
    cache[k]=r; return r
masks=[]
for ln in open(ORBIT):
    mask,size,rep=ln.strip().split('\t');masks.append((int(mask),int(size),rep))
rows=[]
for mask,size,rep in masks:
    grid=[0.25+i*0.025 for i in range(41)] # to 1.25
    vals=[evalp(mask,b) for b in grid]
    br=None
    for a,c in zip(vals,vals[1:]):
        if a['dC']>=0 and c['dC']<=0:
            br=[a['beta'],c['beta']];break
    if br is None:
        best=max(vals,key=lambda r:r['C']); rows.append((mask,size,rep,best['beta'],best['C'],best['dC'],best['skew'],best['binder'],'UNBRACKETED'));continue
    lo,hi=br
    for _ in range(24):
        mid=(lo+hi)/2; r=evalp(mask,mid)
        if r['dC']>0:lo=mid
        else:hi=mid
    r=evalp(mask,(lo+hi)/2)
    rows.append((mask,size,rep,r['beta'],r['C'],r['dC'],r['skew'],r['binder'],'BRACKETED'))
with open(OUT,'w') as f:
    f.write('mask\torbit_size\trep\tbeta_peak\tC_peak\tdC\tskew\tbinder\tstatus\n')
    for r in rows:f.write('\t'.join(map(str,r))+'\n')
print(open(OUT).read())
