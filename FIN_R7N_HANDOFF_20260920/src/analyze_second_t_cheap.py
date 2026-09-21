from pathlib import Path
from fractions import Fraction as F
import json,sys,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); H=ROOT/'inputs/FR223_20260916'; sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
import off_face
IN=ROOT/'checkpoints/R7N-021_t_refine_second_cheap_v1.json'; PREV=ROOT/'checkpoints/R7N-021_t_refine_once_v2.json'; OUT=ROOT/'results/R7N-021_second_t_analysis.json'
def parse(c): return tuple(tuple(map(F,p)) for p in c)
def vol(c):
 v=F(1)
 for a,b in c:v*=b-a
 return v
def l2(c):
 x=np.array([float((a+b)/2) for a,b in c]);p=off_face.p_from_aligned_compact(*x);_,_,_,_,C4=off_face.constants();mu=p@C4;Y=C4-mu;M=Y.T@(p[:,None]*Y);return float(np.linalg.eigvalsh(M)[-2])
d=json.load(open(IN));prev=json.load(open(PREV));h=parse(d['root_hull']);hv=vol(h)
rv=sum((vol(parse(x['cell'])) for x in d['refined_failed']),F(0));prv=sum((vol(parse(x['cell'])) for x in prev['refined_failed']),F(0));sv=sum((vol(parse(x['cell'])) for x in d['refined_safe']),F(0))
vals=[l2(parse(x['cell'])) for x in d['refined_failed']]
out={'task':'R7N-021-second-t-analysis','complete':True,'parent_count':d['parent_residual_count'],'child_count':len(d['refined_safe'])+len(d['refined_failed']),
'safe_count':len(d['refined_safe']),'residual_count':len(d['refined_failed']),'previous_residual_volume_fraction':float(prv/hv),'new_residual_volume_fraction':float(rv/hv),'absolute_volume_reduction_fraction_of_hull':float((prv-rv)/hv),'relative_residual_volume_reduction':float((prv-rv)/prv),'safe_volume_fraction_of_hull':float(sv/hv),'residual_center_lambda2':{'min':min(vals),'median':float(np.median(vals)),'p95':float(np.quantile(vals,.95)),'max':max(vals),'count_over_0_26':sum(z>0.26 for z in vals),'count_over_tau':sum(z>0.268 for z in vals)}}
OUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
