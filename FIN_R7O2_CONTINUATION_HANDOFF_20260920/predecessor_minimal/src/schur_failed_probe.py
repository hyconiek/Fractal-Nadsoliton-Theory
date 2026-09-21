from pathlib import Path
from fractions import Fraction as F
from functools import reduce
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path[:0]=[str(ROOT/'src'),str(ROOT/'inputs/FR223_20260916/src')]
import global_schur_shifted as gs
TAU=F(67,250)
D=json.load(open(ROOT/'checkpoints/R7N-020_trace_e2_compression_v1.json'))
def parsec(x):return [tuple(map(F,p)) for p in x['cell']]
def vol(c):return reduce(lambda a,b:a*b,[hi-lo for lo,hi in c],F(1))
arr=sorted(D['compression_failed'],key=lambda x:vol(parsec(x)),reverse=True)
n=int(sys.argv[1]) if len(sys.argv)>1 else 20
out=[];st=time.time()
for x in arr[:n]:
 B=tuple(tuple(map(F,p)) for p in x['local_box'])
 try:r=gs.raw_shifted_box(B,TAU)
 except Exception as ex:r={'status':'ERROR','reason':type(ex).__name__+':'+str(ex)}
 out.append({'path':x['path'],'status':r['status'],'reason':r.get('reason'),
             'P_hi':float(r['P'].hi) if 'P' in r else None,'P1_lo':float(r['P1'].lo) if 'P1' in r else None,'c2_lo':float(r['c2'].lo) if 'c2' in r else None})
res={'task':'R7N-020-schur-probe','n':n,'pass':sum(x['status']=='INTERVAL_CERTIFIED' for x in out),'error':sum(x['status']=='ERROR' for x in out),'elapsed_seconds':time.time()-st,'results':out}
json.dump(res,open(ROOT/'results/R7N-020_schur_probe.json','w'),indent=2);print(json.dumps({k:v for k,v in res.items() if k!='results'},indent=2))
