from pathlib import Path
from fractions import Fraction as F
import sys, json, hashlib, re
import numpy as np
ROOT=Path('/mnt/data/r7n_repo_root'); CAMP=Path('/mnt/data/fin_rank7_next_campaign')
sys.path.insert(0,str(ROOT))
from fin_rank7_intake_review import scientific_rechecks as s

def compat_I(x):
    if isinstance(x,np.generic): x=x.item()
    if isinstance(x,str) and x.startswith('np.float64(') and x.endswith(')'):
        x=x[len('np.float64('):-1]
    f=F(str(x)) if not isinstance(x,F) else x
    return s.iv.mpf(f.numerator)/f.denominator
s.I=compat_I
s.phases()
new=json.loads((ROOT/'fin_rank7_intake_review/phase_recertification.json').read_text())
old=json.loads((Path('/mnt/data/fin_rank7_intake_review_unpacked/fin_rank7_intake_review/phase_recertification.json')).read_text())
# Normalize only numpy2 spelling inside diagnostic preconditioner strings.
def norm(obj):
    if isinstance(obj,dict): return {k:norm(v) for k,v in obj.items()}
    if isinstance(obj,list): return [norm(v) for v in obj]
    if isinstance(obj,str):
        m=re.fullmatch(r'np\.float64\(([-+0-9.eE]+)\)',obj)
        if m:return m.group(1)
    return obj
nn,nold=norm(new),norm(old)
summary={
 'quartic_count':new['quartic']['count'],'full_count':new['full']['count'],
 'quartic_index_histogram':{},'full_index_histogram':{},
 'normalized_exact_match_to_accepted_20260919':nn==nold,
 'raw_sha256':hashlib.sha256(json.dumps(new,sort_keys=True).encode()).hexdigest()
}
for kind in ['quartic','full']:
    h={}
    for r in new[kind]['roots']:h[str(r['negative_index'])]=h.get(str(r['negative_index']),0)+1
    summary[kind+'_index_histogram']=h
(CAMP/'results/R7N-035_041_phase_local_recertification.json').write_text(json.dumps(new,indent=2)+'\n')
(CAMP/'results/R7N-035_041_phase_replay_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
