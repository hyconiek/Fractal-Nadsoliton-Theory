from pathlib import Path
from fractions import Fraction as F
import sys, json, time
import numpy as np
ROOT=Path('/mnt/data/r7n_repo_root')
sys.path.insert(0,str(ROOT))
from fin_rank7_intake_review import scientific_rechecks as s

def compat_I(x):
    if isinstance(x,np.generic): x=x.item()
    if isinstance(x,str) and x.startswith('np.float64(') and x.endswith(')'):
        x=x[len('np.float64('):-1]
    f=F(str(x)) if not isinstance(x,F) else x
    return s.iv.mpf(f.numerator)/f.denominator
s.I=compat_I
src=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'
for kind,path in [('quartic','results/R7P-089_quartic_roots.json'),('full','certificates/R7P-092_full_phase_roots.json')]:
    data=json.loads((src/path).read_text())['roots']
    t=time.monotonic(); out=s.phase_root(data[0],kind)
    print(kind,'one_root_seconds',time.monotonic()-t,'index',out['negative_index'],'contraction',out['contraction_upper'])
