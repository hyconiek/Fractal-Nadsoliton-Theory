from pathlib import Path
from fractions import Fraction as F
import sys,json
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import boundary_ising as bi
L=bi.strict_intervals(); l3,l4,l5=L[3],L[4],L[5]
sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
tau=F(67,250); gmax=F(250,67); glocal=F('3.71834489812')
out={'task':'R7N-017','tau0':str(tau),'g_endpoint':str(gmax),'sigma_interval':[str(sigma.lo),str(sigma.hi)],
     'sigma_lt_tau0':bool(sigma.hi<tau),'sigma_separation_lower':str(tau-sigma.hi),
     'local_event_gain_decimal_exact':str(glocal),'local_event_below_gain_endpoint':bool(glocal<gmax),
     'gain_margin_exact':str(gmax-glocal),
     'consequence':'If Target P is proved globally, H4=I/g-M4 has at most one strictly negative eigenvalue for 0<g<=250/67. At the endpoint a non-strict ceiling does not imply nonsingularity.',
     'nonconclusions':['not Target S','not X7','not a stationary census','not global minimizer uniqueness']}
json.dump(out,open(ROOT/'results/R7N-017_target_p_implication.json','w'),indent=2)
print(json.dumps(out,indent=2))
