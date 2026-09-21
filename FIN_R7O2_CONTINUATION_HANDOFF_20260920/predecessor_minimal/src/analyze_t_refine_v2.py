from pathlib import Path
from fractions import Fraction as F
import json,sys,math,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); H=ROOT/'inputs/FR223_20260916'; sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
import off_face
IN=ROOT/'checkpoints/R7N-021_t_refine_once_v2.json'; PREV=ROOT/'checkpoints/R7N-021_refine_once_v1.json'
OUT=ROOT/'results/R7N-021_t_refine_analysis.json'

def parse(c): return tuple(tuple(map(F,p)) for p in c)
def vol(c):
 v=F(1)
 for a,b in c:v*=b-a
 return v

def l2_center(c):
 x=np.array([float((a+b)/2) for a,b in c])
 p=off_face.p_from_aligned_compact(*x)
 _,_,_,_,C4=off_face.constants(); mu=p@C4; Y=C4-mu; M=Y.T@(p[:,None]*Y)
 w=np.linalg.eigvalsh(M)
 return float(w[-2]), [float(z) for z in w]

d=json.load(open(IN)); assert d['complete'] and not d['pending_parents']
hull=parse(d['root_hull']); hv=vol(hull)
res=[parse(x['cell']) for x in d['refined_failed']]; rv=sum((vol(c) for c in res),F(0))
safe=[parse(x['cell']) for x in d['refined_safe']]; sv=sum((vol(c) for c in safe),F(0))
prev=json.load(open(PREV)); prv=sum((vol(parse(x['cell'])) for x in prev['refined_failed']),F(0))
l2=[]
for c in res:
 z,e=l2_center(c); l2.append(z)
width_stats=[]
for ax,name in enumerate(['r','s','t','y']):
 vals=[float(b-a) for c in res for a,b in [c[ax]]]
 ratios=[float(b/a) for c in res for a,b in [c[ax]]]
 width_stats.append({'axis':name,'min_width':min(vals),'median_width':float(np.median(vals)),'max_width':max(vals),
                     'min_ratio':min(ratios),'median_ratio':float(np.median(ratios)),'max_ratio':max(ratios)})
reasons={}
for x in d['refined_safe']: reasons[x['reason']]=reasons.get(x['reason'],0)+1
out={'task':'R7N-021-t-refine-analysis','complete':True,'threshold':'67/250',
     'parent_count':d['parent_residual_count'],'child_count':len(d['refined_safe'])+len(d['refined_failed']),
     'safe_child_count':len(d['refined_safe']),'safe_reasons':reasons,'residual_child_count':len(res),
     'previous_residual_volume_fraction':float(prv/hv),'new_residual_volume_fraction':float(rv/hv),
     'absolute_volume_reduction_fraction_of_hull':float((prv-rv)/hv),
     'relative_residual_volume_reduction':float((prv-rv)/prv),
     't_refined_safe_volume_fraction_of_hull':float(sv/hv),
     'residual_center_lambda2':{'min':min(l2),'median':float(np.median(l2)),'p95':float(np.quantile(l2,.95)),'max':max(l2),'count_over_0_26':sum(x>0.26 for x in l2),'count_over_tau':sum(x>0.268 for x in l2)},
     'residual_widths':width_stats}
OUT.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
