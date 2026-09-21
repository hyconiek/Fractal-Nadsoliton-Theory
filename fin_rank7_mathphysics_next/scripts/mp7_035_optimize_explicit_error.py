#!/usr/bin/env python3
import json, math, os, pathlib
import numpy as np
WORK=pathlib.Path(os.environ.get('MP7_WORK_ROOT',pathlib.Path(__file__).resolve().parents[1]))
R7P=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
scan=json.load(open(WORK/'results/MP7-035_cap_lipschitz_scan.json'))
w=json.load(open(WORK/'results/MP7-034_local_phase_weights.json'))
c=json.load(open(R7P/'certificates/R7P-026_equal_energy_event.json'))
glo=float(c['root_box'][4][0])
sqrt_detH={
 'localized':math.sqrt(w['det_G7_localized'][1]/glo**7),
 'uniform':math.sqrt(w['det_G7_uniform'][1]/glo**7),
}

def relerr(cap,sdh,N,cscale):
 m=cap['lipschitz_lambda_min_lower']; B=cap['T3_frob_upper']; r=cap['radius']; d=7
 rho=cscale*math.sqrt(math.log(N)/N)
 if rho>r:return None
 delta=B*cscale**3*math.log(N)**1.5/(6*math.sqrt(N))
 tail=2*d*math.exp(-N*m*rho*rho/(2*d))
 qg=min(1.0,tail)
 C=sdh/(m**(d/2))
 qa=C*tail

 if delta>700:
  return 1e300, {'rho':rho,'delta':delta,'qg':qg,'qa':qa,'ratio_lower':0.0,'ratio_upper':1e300}
 lower=max(0.0,math.exp(-delta)*(1-qg)); upper=math.exp(delta)+qa
 return max(1-lower,upper-1), {'rho':rho,'delta':delta,'qg':qg,'qa':qa,'ratio_lower':lower,'ratio_upper':upper}

def best_for(name,target):
 best=None
 for cap in scan[name]:
  if cap['lipschitz_lambda_min_lower']<=0:continue
  m=cap['lipschitz_lambda_min_lower']
  for expo in np.linspace(0.5,5.0,91):
   cs=math.sqrt(14*expo/m)
   # log scan then binary
   found=None
   for logN in np.linspace(4,25,421):
    N=10**logN; z=relerr(cap,sqrt_detH[name],N,cs)
    if z and z[0]<=target:found=N;break
   if found is None:continue
   a=found/10**0.05;b=found
   for _ in range(50):
    mid=math.sqrt(a*b);z=relerr(cap,sqrt_detH[name],mid,cs)
    if z and z[0]<=target:b=mid
    else:a=mid
   err,detail=relerr(cap,sqrt_detH[name],b,cs)
   cand={'N0':b,'cscale':cs,'tail_exponent':expo,'cap':cap,'detail':detail,'error':err}
   if best is None or b<best['N0']:best=cand
 return best
out={'task':'MP7-035 optimized explicit fixed-cap Laplace error','scientific_state':'PROVED_INTERVAL_ASSISTED_LOCAL_ERROR','targets':{}}
for t in [0.25,0.10,0.05,0.01]:
 bl=best_for('localized',t); bu=best_for('uniform',t); N=max(bl['N0'],bu['N0'])
 # reevaluate both at common N
 el,dl=relerr(bl['cap'],sqrt_detH['localized'],N,bl['cscale']); eu,du=relerr(bu['cap'],sqrt_detH['uniform'],N,bu['cscale'])
 out['targets'][str(t)]={'N0_both_caps':N,'localized':{**bl,'at_common_N_error':el,'at_common_N_detail':dl},'uniform':{**bu,'at_common_N_error':eu,'at_common_N_detail':du}}
out['scope']='fixed caps at the certified equal-energy gain only; no quantitative global-complement error and no uniform moving-minimum g_eq+c/N error'
print(json.dumps(out,indent=2))
