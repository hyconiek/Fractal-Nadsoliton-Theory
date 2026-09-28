#!/usr/bin/env python3
from pathlib import Path
import sys, json, math, csv, hashlib
import numpy as np
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'PHYS-002'))
import safe_seed as s

_,_,laplam,X,A=s.build_rank7(); eig=np.fft.fft(A[0]).real; tr=float(np.trace(A))

def dump(p,obj): p.write_text(json.dumps(obj,indent=2,sort_keys=True),encoding='utf-8')
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tv(p,q): return float(.5*np.abs(p-q).sum())
def kl(p,q): return float(np.sum(p*np.log(p/q)))
def chernoff(p,q):
    r=minimize_scalar(lambda z: math.log(float(np.sum((p**z)*(q**(1-z))))),bounds=(0,1),method='bounded')
    return float(-r.fun),float(r.x)

# A7 structural facts
struct={
 'symmetry_max_abs':float(np.max(abs(A-A.T))),
 'centering_max_abs':float(np.max(abs(A.sum(1)))),
 'diag_spread':float(np.ptp(np.diag(A))),
 'rank_tol_1e-10':int(np.linalg.matrix_rank(A,tol=1e-10)),
 'min_eigenvalue':float(np.linalg.eigvalsh(A).min()),
 'trace':tr,
 'fourier_eigenvalues_k0_6':[float(x) for x in eig[:7]],
 'active_ratios_to_k3':{f'L{k}/L3':float(eig[k]/eig[3]) for k in [4,5,6]},
 'construction_laplacian_modes_k0_6':[float(x) for x in laplam],
}

# PHYS-001 enumeration identity tests
p1={}
for N in [1,2]:
    states=s.compositions(N); pi_count=s.count_pi(states,0.731,A,theta=.17)
    # enumerate labelled microstates, aggregate by counts
    agg={st:0.0 for st in states}; raw=[]
    for sig in __import__('itertools').product(range(s.Q), repeat=N):
        n=np.bincount(sig,minlength=s.Q).astype(float)
        Hexpo=(.731/(2*N))*float(n@A@n)+.17*n[0]
        raw.append((sig,Hexpo,tuple(map(int,n))))
    z=sum(math.exp(x[1]) for x in raw)
    for sig,lw,st in raw: agg[st]+=math.exp(lw)/z
    err=max(abs(agg[st]-pi_count[i]) for i,st in enumerate(states))
    # conditional check for one context when possible
    p1[str(N)]={'count_states':len(states),'labelled_states':s.Q**N,'max_aggregated_probability_error':float(err)}
p1['N2_conditional_row0_max_error']=float(np.max(abs(s.pair_difference_distribution(A,.731)-
    np.exp((.731/2)*A[0]-np.max((.731/2)*A[0]))/np.exp((.731/2)*A[0]-np.max((.731/2)*A[0])).sum())))
p1['structural']=struct
dump(ROOT/'PHYS-001'/'enumeration_tests.json',p1)

# PHYS-002 diagnostics
repro={'baseline_commit':'97c4231f33632800fd817fe2294555bb8bcb041f','A7':struct,'generators_N2':{},'sector_g0_N2':{},'rho3':{},'source_gap':{'historical_FIN_son_files_tracked_at_baseline':False}}
for kin in ['heat_bath','metropolis','barker']:
    st=s.compositions(2); pi=s.count_pi(st,s.G_FROZEN,A); q=s.count_generator(st,s.G_FROZEN,A,kin)
    chk=s.generator_checks(st,q,pi); rev,secs=s.sector_spectra(st,q,pi)
    chk['reversible_similarity_symmetry_defect']=rev; chk['max_sector_eigenpair_residual']=max(v['max_eigenpair_residual'] for v in secs.values())
    repro['generators_N2'][kin]=chk
st=s.compositions(2); pi=s.count_pi(st,0,A); q=s.count_generator(st,0,A); rev,secs=s.sector_spectra(st,q,pi)
repro['sector_g0_N2']={'all_sector_dimensions':{k:v['dimension'] for k,v in secs.items()},'max_eigenpair_residual':max(v['max_eigenpair_residual'] for v in secs.values()),'all_12_sectors_present':all(v['dimension']>0 for v in secs.values())}
st=s.compositions(3); pi=s.count_pi(st,s.G_FROZEN,A); q=s.count_generator(st,s.G_FROZEN,A); rev,secs=s.sector_spectra(st,q,pi)
rho=-max(secs['4']['eigenvalues']); repro['rho3']={'computed':rho,'fixture':0.13143978619564534,'abs_error':abs(rho-0.13143978619564534),'sector_k':4,'eigenpair_residual':secs['4']['max_eigenpair_residual']}
repro['G_FROZEN']={'code_value':s.G_FROZEN,'R09_g_bal':5.145228719489144,'difference':s.G_FROZEN-5.145228719489144,'role':'operational barrier-balance value, not physical constant'}
dump(ROOT/'PHYS-002'/'REPRODUCTION.json',repro)

# PHYS-003 static spectral fingerprints
fp={'formula':{'S0':1.0,'slope':'(N-1)*Lambda_k/(12*N)'},'A7':struct,'smallN':{},'finite_g_remainder':{}}
for N in [1,2,3]:
    st=s.compositions(N)
    fp['smallN'][str(N)]={}
    for g in [0.0,0.05,0.1,0.2]:
        pi=s.count_pi(st,g,A)
        fp['smallN'][str(N)][str(g)]={str(k):s.static_S(st,pi,N,k) for k in range(1,7)}
    if N>1:
        rem={}
        for g in [.05,.1,.2]:
            vals=fp['smallN'][str(N)][str(g)]
            rr={}
            for k in range(1,7):
                slope=(N-1)*eig[k]/(12*N)
                rr[str(k)]=float(vals[str(k)]-(1+slope*g))
            rem[str(g)]={'per_k':rr,'max_abs':max(abs(v) for v in rr.values())}
        fp['finite_g_remainder'][str(N)]=rem
fp['N2_pair_g_0p2']=s.pair_difference_distribution(A,.2).tolist()
fp['minimal_readout']='microscopic 12-label pair difference (or equivalent Fourier histogram), not basin label'
dump(ROOT/'PHYS-003'/'FINGERPRINT.json',fp)

# PHYS-004 models and static/dynamic contrasts
mods=s.build_countermodels(A)
model_json={}
for name,M in mods.items():
    ev=np.fft.fft(M[0]).real
    model_json[name]={'trace':float(np.trace(M)),'rank':int(np.linalg.matrix_rank(M,tol=1e-10)),'fourier_k0_6':[float(x) for x in ev[:7]],'min_eigenvalue':float(np.linalg.eigvalsh(M).min())}
dump(ROOT/'PHYS-004'/'COUNTERMODELS.json',model_json)
rows=[]
for g in [.2,1.0,2.0,3.0]:
    p=s.pair_difference_distribution(A,g)
    for name,M in mods.items():
        if name=='FIN_A7': continue
        q=s.pair_difference_distribution(M,g); C,cs=chernoff(p,q)
        rows.append({'kind':'static_pair','g':g,'model':name,'TV':tv(p,q),'KL_FIN_to_alt':kl(p,q),'Chernoff':C,'chernoff_s':cs,'M_bound_5pct_equal_prior':math.ceil(math.log(10)/C) if C>0 else ''})
# dynamics for N=2 at G; normalize clocks by k4 rate at g=0
for kin in ['heat_bath','metropolis','barker']:
    st=s.compositions(2); pi0=s.count_pi(st,0,A); q0=s.count_generator(st,0,A,kin); _,sec0=s.sector_spectra(st,q0,pi0); r0=-max(sec0['4']['eigenvalues'])
    pi=s.count_pi(st,s.G_FROZEN,A); q=s.count_generator(st,s.G_FROZEN,A,kin); _,sec=s.sector_spectra(st,q,pi)
    for k in range(1,7):
        rate=-max(sec[str(k)]['eigenvalues'])/r0
        rows.append({'kind':'dynamic_sector','g':s.G_FROZEN,'model':kin,'sector_k':k,'normalized_rate':rate,'time_normalization':'divide by k4 rate at g=0'})
with open(ROOT/'PHYS-004'/'CONTRAST_MATRIX.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=sorted(set().union(*[r.keys() for r in rows]))); w.writeheader(); w.writerows(rows)
dump(ROOT/'PHYS-004'/'contrast_summary.json',{'rows':rows})

# PHYS-007 prereg predictions at g=3 pair-difference test
GTEST=3.0; p=s.pair_difference_distribution(A,GTEST)
contr={}
for name,M in mods.items():
    if name=='FIN_A7': continue
    q=s.pair_difference_distribution(M,GTEST); C,cs=chernoff(p,q)
    contr[name]={'TV':tv(p,q),'KL':kl(p,q),'Chernoff':C,'M_equal_prior_bound_5pct':math.ceil(math.log(10)/C) if C>0 else None}
# calibration sensitivity +/-0.2% g
sens=max(tv(p,s.pair_difference_distribution(A,GTEST*(1+d))) for d in [-.002,.002])
pr={'test':'N=2 conditional pair-difference histogram with anchor label i=0','N':2,'g_test':GTEST,'g_zero_control':0.0,'calibration_g_relative_tolerance':0.002,'max_TV_shift_from_g_tolerance':sens,'primary_countermodels':['FULL_POTTS_TRACE','FLAT_P7_TRACE','PERT_34_10P','PERT_35_10P'],'sensitivity_only_countermodels':['PERT_36_5P','LEAK_K12_2P_TRACE'],'minimum_primary_TV':min(contr[n]['TV'] for n in ['FULL_POTTS_TRACE','FLAT_P7_TRACE','PERT_34_10P','PERT_35_10P']),'independent_samples':'independent reset/anchor cycles; trajectory data require ESS correction','burn_in':'not required if one site is clamped and the other is drawn directly from the calibrated conditional; otherwise must be validated separately','no_validation_records_opened':True}
dump(ROOT/'PHYS-007'/'PREREGISTRATION.json',pr)
dump(ROOT/'PHYS-007'/'PREDICTIONS.json',{'FIN_pair_difference':p.tolist(),'countermodel_metrics':contr})
dump(ROOT/'PHYS-007'/'CONTRASTS.json',contr)

print(json.dumps({'A7':struct,'rho3':repro['rho3'],'primary_min_TV_g3':pr['minimum_primary_TV'],'g_uncertainty_TV':sens,'contrasts_g3':contr},indent=2))
