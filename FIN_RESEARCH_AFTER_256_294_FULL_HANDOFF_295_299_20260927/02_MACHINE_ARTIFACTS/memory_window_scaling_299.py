#!/usr/bin/env python3
import json, math, hashlib
from pathlib import Path
import numpy as np

OUT=Path('/mnt/data/fin299')
SOURCES={
 'crt_slow_memory_scaling_N3_N8.json':'adc611c5da2bad0bbe73e7f3b56331cf58a1976d',
 'largeN_memory_N7.json':'a36dd3dfb6258bf82a481bf96694c69910c08af0',
 'largeN_memory_N8.json':'03b15853dd4cb807cd6ce525428f628a0b11c28d',
 'rigorous_memory_tail_bounds_N6.json':'3c5923466d1da889e5e9b1ff602f98adb75cde1b',
 'STIELTJES_MEMORY_STRUCTURE_231.md':'15d3f314c7fe439dee76fa64950cccef3c20e987',
}
rows=[
 dict(N=3,rho=0.13200595666677714,lambda_exact=-0.13143978619564634,M0=0.23725111000188426,M1=0.08822785603151304,mz_relerr=0.004307451248346287),
 dict(N=4,rho=0.07185438302832774,lambda_exact=-0.0716922774819137,M0=0.214893967498199,M1=0.07729026831880151,mz_relerr=0.00226112981910681),
 dict(N=5,rho=0.04022475665704753,lambda_exact=-0.040186755916096495,M0=0.11687202466222309,M1=0.04142515358374989,mz_relerr=0.0009456035971246961),
 dict(N=6,rho=0.02261891215470609,lambda_exact=-0.02261127518711436,M0=0.09063392463179727,M1=0.028050238846990486,mz_relerr=0.0003377504156016344),
 dict(N=7,rho=0.012641911063094438,lambda_exact=-0.012640386432677186,M0=0.04826980176722458,M1=0.014816161830710436,mz_relerr=0.0001206158075445394),
 dict(N=8,rho=0.006989922026846724,lambda_exact=-0.00698967186661201,M0=0.031939601938700726,M1=0.008650007293639801,mz_relerr=3.578998263262173e-05),
]
for r in rows:
    r['tau_mem_mean']=r['M1']/r['M0']
    r['tau_slow']=1/r['rho']
    r['epsilon_moment']=r['rho']*r['tau_mem_mean']
    r['epsilon_moment_exact_clock']=abs(r['lambda_exact'])*r['tau_mem_mean']
    r['Z_onepole']=1/(1+r['M1'])
    r['initial_slip_deficit']=1-r['Z_onepole']
    r['markov_tail_bound_theta_0p05']=min(1.0,r['epsilon_moment']/0.05)
    r['markov_tail_bound_theta_0p10']=min(1.0,r['epsilon_moment']/0.10)
    r['markov_tail_bound_theta_0p20']=min(1.0,r['epsilon_moment']/0.20)

# Coupled hidden spectral gaps in the slow k=4 / Z3 sector. N=6 comes from the proof-grade tail-bound artifact;
# N=7,8 are finite-state numerical eigensolve values from the accepted large-N memory artifacts.
gap_data={6:0.6771882333363985,7:0.484278766376313,8:0.4997526324109206}
for r in rows:
    if r['N'] in gap_data:
        g=gap_data[r['N']]
        r['coupled_hidden_gap']=g
        r['epsilon_gap']=r['rho']/g
        for theta in (0.05,0.10,0.20):
            r[f'spectral_tail_bound_theta_{theta:.2f}']=math.exp(-theta/r['epsilon_gap'])

Ns=np.array([r['N'] for r in rows],float)
eps=np.array([r['epsilon_moment'] for r in rows],float)
mz=np.array([r['mz_relerr'] for r in rows],float)

def exp_fit(y):
    c=np.polyfit(Ns,np.log(y),1)
    beta=-float(c[0]); A=float(np.exp(c[1])); pred=A*np.exp(-beta*Ns)
    r2=1-float(np.sum((np.log(y)-np.log(pred))**2)/np.sum((np.log(y)-np.mean(np.log(y)))**2))
    return {'A':A,'beta_per_N':beta,'factor_per_N':math.exp(-beta),'log_R2':r2}

fits={'epsilon_moment_descriptive_exponential':exp_fit(eps),'mz_relerr_descriptive_exponential':exp_fit(mz)}

# Exact scalar Stieltjes-tail theorem audit:
# f(t)=K(t)/M0 is a probability density, E_f[t]=M1/M0.
# Markov: int_T^inf K / M0 <= (M1/M0)/T.
# If coupled support gamma>=gamma_c: int_T^inf K / M0 <= exp(-gamma_c T).
assert all(rows[i+1]['epsilon_moment'] < rows[i]['epsilon_moment'] for i in range(len(rows)-1))
assert all(rows[i+1]['tau_slow'] > rows[i]['tau_slow'] for i in range(len(rows)-1))
assert all(rows[i+1]['mz_relerr'] < rows[i]['mz_relerr'] for i in range(len(rows)-1))
assert all(rows[i+1]['initial_slip_deficit'] < rows[i]['initial_slip_deficit'] for i in range(len(rows)-1))

ratio_eps=rows[0]['epsilon_moment']/rows[-1]['epsilon_moment']
ratio_slow=rows[-1]['tau_slow']/rows[0]['tau_slow']
ratio_tau_mem=rows[0]['tau_mem_mean']/rows[-1]['tau_mem_mean']

summary={
 'report':299,
 'status':'finite-N scale-separation certificate plus exact Stieltjes tail inequalities; no N->infinity theorem',
 'source_git_blob_shas':SOURCES,
 'theorem':{
   'positive_memory_kernel':'K(t)=int exp(-gamma t) dmu(gamma), dmu>=0',
   'normalized_density':'f(t)=K(t)/M0',
   'mean_memory_time':'tau_mem=M1/M0',
   'moment_tail_bound':'Tail(T)=(1/M0) int_T^inf K(t)dt <= tau_mem/T',
   'slow_fraction_form':'At T=theta/rho, Tail <= epsilon_moment/theta, epsilon_moment=rho*M1/M0',
   'spectral_tail_bound':'If all coupled hidden rates gamma>=gamma_c, Tail(T)<=exp(-gamma_c*T)=exp(-theta/epsilon_gap), epsilon_gap=rho/gamma_c'
 },
 'rows':rows,
 'fits':fits,
 'finite_range_summary':{
   'epsilon_moment_N3':rows[0]['epsilon_moment'],
   'epsilon_moment_N8':rows[-1]['epsilon_moment'],
   'epsilon_drop_factor_N3_to_N8':ratio_eps,
   'slow_time_growth_factor_N3_to_N8':ratio_slow,
   'mean_memory_time_change_factor_N3_to_N8':ratio_tau_mem,
   'N8_moment_tail_bound_after_10pct_slow_time':rows[-1]['markov_tail_bound_theta_0p10'],
   'N8_spectral_tail_bound_after_10pct_slow_time':rows[-1]['spectral_tail_bound_theta_0.10'],
   'N8_spectral_tail_bound_after_20pct_slow_time':rows[-1]['spectral_tail_bound_theta_0.20'],
 },
 'scope_boundary':[
   'N=3..8 only for moment scaling; N=6..8 only for coupled-gap check',
   'descriptive exponential fits are not asymptotic laws or barrier exponents',
   'N=7,8 coupled gaps are finite-state numerical eigensolve values, not interval certificates',
   'the result supports controlled effective Markov closure but does not derive spatial sites, SI time, or a fundamental FIN completeness law',
   'no QW-2191, legacy-to-strict bridge, Standard Model, gravity, or ToE closure follows'
 ]
}
(OUT/'memory_window_scaling_299.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({
 'PASS':True,
 'eps_N3':rows[0]['epsilon_moment'],
 'eps_N8':rows[-1]['epsilon_moment'],
 'drop_factor':ratio_eps,
 'tau_mem_N3':rows[0]['tau_mem_mean'],
 'tau_mem_N8':rows[-1]['tau_mem_mean'],
 'slow_N3':rows[0]['tau_slow'],
 'slow_N8':rows[-1]['tau_slow'],
 'N8_tail_10pct_moment':rows[-1]['markov_tail_bound_theta_0p10'],
 'N8_tail_10pct_spectral':rows[-1]['spectral_tail_bound_theta_0.10'],
 'eps_fit_factor_per_N':fits['epsilon_moment_descriptive_exponential']['factor_per_N'],
 'eps_fit_log_R2':fits['epsilon_moment_descriptive_exponential']['log_R2'],
 'mz_fit_factor_per_N':fits['mz_relerr_descriptive_exponential']['factor_per_N'],
},indent=2))
