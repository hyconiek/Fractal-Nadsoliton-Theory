from pathlib import Path
import json, math
from fractions import Fraction as F
from .face_certificate import cert_data
from .intervals import QI,sqrt_interval
ROOT=Path(__file__).resolve().parents[1]

def build():
 D=cert_data();L=D['L'];t2=QI(1)-D['sigma']/D['a']
 d3=L[3]*(2*t2-QI(1))/6
 rad=(L[5]-L[4])**2+4*L[4]*L[5]*t2
 denom=6*sqrt_interval(rad)
 # negative quotient: construct with interval ops
 mag=L[4]*L[5]*t2/denom
 d45=-mag
 out={'id':'R7P-023-parity-first-order',
      'claim_id':'CLM-PARITY-SPLIT',
      'domain':'double-root parity-mixing perturbation at q=1',
      'quantifiers':'for every strict spectral tuple in the accepted outward intervals',
      'assumptions':['accepted strict spectral intervals','fixed-t degenerate perturbation at the certified extreme-face double root'],
      'proof_type':'exact symbolic first-order formulas + rational interval sign enclosure',
      'inputs':['inputs/fin_handoff_audit/results.json'],
      'conclusion':'both fixed-t first-order eigenvalue shifts are strictly negative; the reoptimized envelope slope remains uncertified',
      'global_pass':False,
      't_star_squared':[str(t2.lo),str(t2.hi)],
      'delta_k3':[str(d3.lo),str(d3.hi)],
      'delta_45':[str(d45.lo),str(d45.hi)],
      'both_strictly_negative':d3.hi<0 and d45.hi<0,
      'reoptimized_slope_status':'NUMERICAL_SEED_ONLY; imported ~ -0.1312828584 per (1-q), not certified here.'}
 (ROOT/'certificates/R7P-023_parity_first_order.json').write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':print(json.dumps(build(),indent=2))
