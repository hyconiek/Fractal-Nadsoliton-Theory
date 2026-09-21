from fractions import Fraction as F
import sys,unittest,numpy as np
from pathlib import Path
ROOT=Path('/mnt/data/r7n_repo_root');sys.path.insert(0,str(ROOT))
from fin_rank7_intake_review import scientific_rechecks as s

def compat_I(x):
    if isinstance(x,np.generic): x=x.item()
    if isinstance(x,str) and x.startswith('np.float64(') and x.endswith(')'): x=x[len('np.float64('):-1]
    f=F(str(x)) if not isinstance(x,F) else x
    return s.iv.mpf(f.numerator)/f.denominator
s.I=compat_I
from fin_rank7_intake_review import test_integration
suite=unittest.defaultTestLoader.loadTestsFromModule(test_integration)
r=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(0 if r.wasSuccessful() else 1)
