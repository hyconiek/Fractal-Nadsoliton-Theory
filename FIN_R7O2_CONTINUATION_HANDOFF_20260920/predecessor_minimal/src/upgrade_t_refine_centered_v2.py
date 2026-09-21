from pathlib import Path
from fractions import Fraction as F
import sys,json,time,shutil
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); sys.path.insert(0,str(ROOT/'src'))
import compression_centered_moment as cm
SRC=ROOT/'checkpoints/R7N-021_t_refine_once_v1.json'
OUT=ROOT/'checkpoints/R7N-021_t_refine_once_v2.json'

def parse(cell): return tuple(tuple(map(F,p)) for p in cell)
def main():
 d=json.load(open(SRC)); st=time.time(); remain=[]; gained=[]
 for x in d['refined_failed']:
    r=cm.certify(parse(x['cell']))
    if r['ok']:
        gained.append({**x,'reason':'SAFE_BY_CENTERED_MOMENT','centered_reason':r['reason'],
                       'basis_den':r['basis_den'],'basis_num':r['basis_num'],'rank_rows':r['rank_rows'],'rank_det':r['rank_det'],
                       'center_c':r['center_c'],'d1':r['d1'],'d2':r['d2'],'d3':r['d3'],'gersh_lower':r['gersh_lower']})
    else: remain.append(x)
 d['task']='R7N-021-t-refine-once-v2'
 d['policy']='one exact binary split in t; classify trace->e2->centered-moment->cell-adaptive compression; v1 completed prefix upgraded by centered-moment'
 d['refined_safe'].extend(gained); d['refined_failed']=remain
 d['stats']['SAFE_BY_CENTERED_MOMENT']=len(gained)
 d['stats']['UNRESOLVED']=len(remain)
 d['upgrade_from_v1']={'checked_previous_failures':len(gained)+len(remain),'new_centered_safe':len(gained),'seconds':time.time()-st}
 OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps(d['upgrade_from_v1']|{'current_safe_total':len(d['refined_safe']),'current_failed':len(d['refined_failed']),'pending':len(d['pending_parents'])},indent=2))
if __name__=='__main__': main()
