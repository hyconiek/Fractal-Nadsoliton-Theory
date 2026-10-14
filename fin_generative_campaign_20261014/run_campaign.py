#!/usr/bin/env python3
"""Independent bounded reconstruction of RG-C01..C05 from the published FIN formulas.
Requires numpy scipy mpmath. Does not import the original checkout; see provenance.
"""
import os, json, math, hashlib, itertools, csv, time, sys
from pathlib import Path
import numpy as np
from scipy.special import logsumexp, softmax
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm
ROOT=Path(__file__).resolve().parent
Q=12; G=5.145228719489142
BASE={'q':Q,'gain':G,'parent_static_M':4,'parent_dynamic_M':3,'rate_attempt_per_copy':1,'theta':0,'kinetics':'leave_one_out_heat_bath','kernel_parameters':{'omega':'0.18575','phi':'0.1625','eta':'1.8'},'parent_normalization_rule':'fixed M even after deletion','observation':'first retained labelled spins, ordered'}
def write_json(p, obj):
    pp=ROOT/p; pp.parent.mkdir(parents=True,exist_ok=True);pp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True,default=lambda x: x.item() if hasattr(x,'item') else str(x))+'\n')
def write_md(p, s):
    pp=ROOT/p;pp.parent.mkdir(parents=True,exist_ok=True);pp.write_text(s+'\n')
def hashfile(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def build_rank7():
    j=np.arange(Q);W=np.zeros((Q,Q))
    for i in range(Q):
        for k in range(Q):
            if i!=k:
                d=min(abs(i-k),Q-abs(i-k))
                W[i,k]=math.cos(.18575*d+.1625)/(1+d**1.8)
    L=np.diag(W.sum(axis=1))-W; eig=np.fft.fft(L[0]).real[:7]
    cols=[]
    for k in (3,4,5):
        cols.extend([np.sqrt(eig[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(eig[k]/6)*np.sin(2*np.pi*k*j/Q)])
    cols.append(np.sqrt(eig[6]/12)*(-1.)**j)
    X=np.column_stack(cols);A=X@X.T
    return W,eig,X,A
def models_from_A(A):
    trace=float(np.trace(A));P=np.eye(Q)-np.ones((Q,Q))/Q
    lv=np.zeros(Q);lv[[3,4,5,6,7,8,9]]=trace/7
    row=np.fft.ifft(lv).real
    flat=np.array([np.roll(row,i) for i in range(Q)])
    return {'FIN_A7':A,'FULL_POTTS_TRACE':P*(trace/11),'FLAT_P7_TRACE':flat}
def labelled_log(s, M, g, A):
    n=np.bincount(s,minlength=Q);return (g/(2*M))*n@A@n
def marginal_F3(A,g):
    state=np.indices((Q,Q,Q)).reshape(3,-1).T
    counts=np.zeros((len(state),Q),int)
    for k in range(3): counts[np.arange(len(state)),state[:,k]]+=1
    pair=(g/8)*np.einsum('bi,ij,bj->b',counts,A,counts)
    corr=logsumexp((g/4)*(counts@A),axis=1)
    return (pair+corr).reshape((Q,Q,Q))
def three_decompose(F):
    F0=float(F.mean()); Fi=[]
    for i in range(3):
        axes=tuple(j for j in range(3) if j!=i)
        z=F.mean(axis=axes)-F0
        shape=[1,1,1];shape[i]=Q;Fi.append(z.reshape(shape))
    Fij=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        k=3-i-j
        zz=F.mean(axis=k,keepdims=True)-F0-Fi[i]-Fi[j]
        Fij.append(zz)
    F2=F0+sum(Fi)+sum(Fij)
    F3=F-F2
    return F2,F3,{'constant':F0,'one_body_inf_norm':float(max(np.max(np.abs(z)) for z in Fi)), 'two_body_inf_norm':float(max(np.max(np.abs(z)) for z in Fij)),'three_body_inf_norm':float(np.max(np.abs(F3))),'components_center_max_abs':float(max(max(abs(F3.mean(axis=i)).ravel()) for i in range(3))),'reconstruction_max_abs':float(np.max(np.abs(F-(F2+F3))))}
def prob(F):return np.exp(F-logsumexp(F))
def tv(p,q):return .5*float(np.sum(np.abs(p-q)))
def osc(r):return float(np.max(r)-np.min(r))
def alt_contrast(F):
    out=0.
    for b in itertools.product((0,1),repeat=3):
        v=(0 if b[0]==0 else 1,0 if b[1]==0 else 2,0 if b[2]==0 else 4)
        out+=(-1)**sum(b)*float(F[v])
    return out
def sequential_direct_2(A,g):
    out1=np.zeros((Q,Q));out2=np.zeros((Q,Q)); M=4
    for a in range(Q):
        for b in range(Q):
            n=np.bincount([a,b],minlength=Q)
            u=[]
            for c in range(Q):
                nc=n.copy();nc[c]+=1
                u.append((g/(2*M))*nc@A@nc+logsumexp((g/M)*(A@nc)))
            out1[a,b]=logsumexp(u)
            v=[labelled_log((a,b,c,d),4,g,A) for c in range(Q) for d in range(Q)]
            out2[a,b]=logsumexp(v)
    return out1,out2

def parent_N3(A,g):
    states=np.indices((Q,Q,Q)).reshape(3,-1).T
    counts=np.zeros((len(states),Q),int)
    for k in range(3): counts[np.arange(len(states)),states[:,k]]+=1
    energy=(g/6)*np.einsum('bi,ij,bj->b',counts,A,counts)
    pi=prob(energy)
    rr=[];cc=[];dd=[];diag=np.zeros(len(states));idx=lambda s:int((s[0]*Q+s[1])*Q+s[2])
    for ix,s in enumerate(states):
        for which in range(3):
            other=counts[ix].copy();other[s[which]]-=1
            rates=softmax((g/3)*(A@other))
            for j in range(Q):
                if j==s[which]:continue
                s2=s.copy();s2[which]=j
                jj=idx(s2);rate=float(rates[j]);rr.append(ix);cc.append(jj);dd.append(rate);diag[ix]-=rate
    rr.extend(range(len(states)));cc.extend(range(len(states)));dd.extend(diag)
    L=csr_matrix((dd,(rr,cc)),shape=(len(states),)*2)
    yi=states[:,0]*Q+states[:,1]
    # equilibrium projection, conditional weights within every observed block
    pY=np.bincount(yi,weights=pi,minlength=Q*Q)
    conditional=pi/pY[yi]
    coarse=np.zeros((Q*Q,Q*Q)); d= L.tocoo()
    for i,j,r in zip(d.row,d.col,d.data):
        if yi[i]!=yi[j]:coarse[yi[i],yi[j]]+=conditional[i]*r
    np.fill_diagonal(coarse,-coarse.sum(axis=1))
    # different declared generator: heat-bath sampler with exact observed marginal stationary law
    gibbs=np.zeros_like(coarse)
    log_marg=np.log(pY).reshape((Q,Q))
    for i in range(Q):
        for j in range(Q):
            ix=i*Q+j
            for i2 in range(Q):
                if i2!=i:gibbs[ix,i2*Q+j]+=softmax(log_marg[:,j])[i2]
            for j2 in range(Q):
                if j2!=j:gibbs[ix,i*Q+j2]+=softmax(log_marg[i,:])[j2]
    np.fill_diagonal(gibbs,-gibbs.sum(axis=1))
    return L,pi,yi,pY,conditional,coarse,gibbs

def dynamics(A,g, full=False):
    L,pi,yi,pY,cond,B,D=parent_N3(A,g)
    row=float(np.max(np.abs(np.asarray(L.sum(axis=1)))))
    off=L.copy().tolil();off.setdiag(0);off=off.tocsr();off.eliminate_zeros()
    neg=float(off.data.min()) if len(off.data) else 0
    station=float(np.max(np.abs(pi@L)))
    flux= L.multiply(pi[:,None]); db=(flux-flux.T).tocoo();balance=float(np.max(np.abs(db.data))) if len(db.data) else 0
    S=np.sqrt(pi); symmetry=L.multiply(S[:,None]).multiply(1/S[None,:]);symdiff=(symmetry-symmetry.T).tocoo();selfadj=float(np.max(np.abs(symdiff.data))) if len(symdiff.data) else 0
    rates0=softmax((g/3)*(A[:,0]+A[:,0]));rates1=softmax((g/3)*(A[:,0]+A[:,1]))
    valid_visible_targets=np.arange(1,Q) # current visible label is zero; j=0 is no CTMC jump
    jw=int(valid_visible_targets[np.argmax(np.abs(rates0[valid_visible_targets]-rates1[valid_visible_targets]))])
    defect=float(abs(rates0[jw]-rates1[jw]))
    y0=0; ids=np.flatnonzero(yi==y0); p0=np.zeros(len(pi));p0[ids]=cond[ids]
    pred0=np.zeros(len(pY));pred0[y0]=1
    times=[.15,.5,1.,1.5] if full else [.5,1.]
    measures=[]
    for t in times:
        exact=np.bincount(yi,weights=expm_multiply(L.T*t,p0),minlength=Q*Q)
        av=expm_multiply(B.T*t,pred0);gf=expm_multiply(D.T*t,pred0)
        measures.append({'t':t,'TV_exact_vs_equilibrium_generator':tv(exact,av),'TV_exact_vs_gibbs_sampler':tv(exact,gf),'TV_two_autonomous_generators':tv(av,gf), 'event00_exact':float(exact[0]),'event00_markov':float(av[0])})
    # Explicit hidden nonequilibrium prep, fixed same y0
    p00=np.zeros(len(pi));p00[0]=1
    p01=np.zeros(len(pi));p01[1]=1
    prep=[]
    for t in [.15,.5]:
        a=np.bincount(yi,weights=expm_multiply(L.T*t,p00),minlength=Q*Q)
        b=np.bincount(yi,weights=expm_multiply(L.T*t,p01),minlength=Q*Q)
        prep.append({'t':t,'TV_hidden_0_vs_hidden_1':tv(a,b),'TV_hidden_0_vs_stationary_markov':tv(a,expm_multiply(B.T*t,pred0))})
    # Quadratic reversible memory channel f=1_{Y=(0,0)}.
    f=(yi==0).astype(float);Lf=L@f
    Plf=np.bincount(yi,weights=pi*Lf,minlength=Q*Q)/pY
    v=Lf-Plf[yi]
    M0=float(np.dot(pi,v*v)); M1=float(np.dot(pi*v,L@v));
    gap=float(np.max(np.abs(np.bincount(yi,weights=pi*v,minlength=Q*Q))))
    out={'generator_checks':{'row_sum_abs':row,'min_off_diagonal':neg,'stationarity_abs':station,'detailed_balance_abs':balance,'reversible_matrix_symmetry_abs':selfadj}, 'model_kinetics':'per-copy CTMC leave-one-out heat_bath at attempt rate 1','strong_lumpability_witness':{'observed_state':[0,0],'hidden_microstates':[[0,0,0],[0,0,1]],'maximum_single_visible_refresh_rate_difference':defect,'j_witness':jw,'rates_hidden0':float(rates0[jw]),'rates_hidden1':float(rates1[jw])},'projection':{'equilibrium_distribution_difference_av_vs_marg':float(np.max(np.abs(pY@B))), 'equilibrium_distribution_difference_gibbs_vs_marg':float(np.max(np.abs(pY@D))), 'max_generator_row_difference':float(np.max(np.sum(np.abs(B-D),axis=1)))},'finite_window':measures,'hidden_preparation_dependence':prep,'memory_quadratic_channel':{'observable':'1{first_two_labels=(0,0)}','M0_coupled':M0,'M1_derivative':M1,'projection_orthogonality_error':gap, 'positivity_scope':'equilibrium reversible quadratic scalar channel'}}
    if full:
        # QLP = v; QLQ generator on hidden fluctuations.
        R=np.arange(Q*Q)
        rr=[];cc=[];dat=[]
        for y in R:
            ix=np.flatnonzero(yi==y)
            for a in ix:
                rr.append(a);cc.append(a);dat.append(1.)
                for b in ix:rr.append(a);cc.append(b);dat.append(-cond[b])
        QP=csr_matrix((dat,(rr,cc)),shape=L.shape)
        H=(QP@L@QP).tocsr()
        out['memory_quadratic_channel']['K_t']={str(t):float(np.dot(pi*v,expm_multiply(H*t,v))) for t in [0,.2,.5,1.]}
        out['memory_quadratic_channel']['memory_initial_slip_nonstationary_prep']=float((L@f)[0]-Plf[0])
    return out

def safe_float(obj):return json.loads(json.dumps(obj,default=lambda x: x.item() if hasattr(x,'item') else str(x)))
def task_docs(task,title,data,report,limits):
    p=task+'/'
    write_json(p+'RESULTS.json',data)
    write_md(p+'REPORT.md','# '+task+' — '+title+'\n\n'+report)
    write_md(p+'NONCONCLUSIONS.md','# Limits / nonconclusions\n\n'+limits)
    write_md(p+'REPLAY.md', '# Reproduction\n\nFrom the top-level campaign directory run `python run_campaign.py`. Requires Python 3, numpy, scipy and mpmath (for optional interval source check). `python independent_verify.py` performs an independently written consistency check. No external datasets or hardware are used. See SOURCE_ACCESS.md for the missing repository snapshot.\n')

def main():
    started=time.time(); W,eig,X,A=build_rank7();models=models_from_A(A)
    current={'execution_status':'EXECUTED_WITH_SOURCE_LIMITATION','scientific_verdict':'FINITE_MODEL_CONTRACT_ONLY','model_id':'FIN_A7_safe_seed_reconstruction_20261014','assumptions':BASE,'supplied_objects':['kernel tuple','rank7 Fourier construction','frozen gain','single-copy heat bath CTMC'], 'derived_objects':['A7 numeric reconstruction','matrix trace and symmetry'], 'fitted_objects':[],'arithmetic_status':'binary64 numerical; separate interval coefficient verification','error_budget':{'numerical_source_precision':'see interval_source_check.json'},'discrepancies':['full Git checkout and commit timeline unavailable in execution container'],'dependencies':['public safe_seed.py', 'AGENTS.md', 'attached plan'],'next_action':'perform finite static elimination and dynamic checks'}
    write_json('RG-C01/MODEL_CONTRACT.json',{'contract':BASE,'matrix_A7':A.tolist(),'lambda_k_first7':eig.tolist(),'source_label':'independent transcription of source formula; NOT file-identical repository safe_seed.py','comparison_models':list(models),'full_source_commit':'NOT_VERIFIED','declared_baseline_in_plan':'c24773ec'})
    write_json('RG-C01/INPUT_REGISTRY.json', {'actual_local_input':['provided markdown plan'],'verified_remote_via_web':['AGENTS.md','PHYS-002/safe_seed.py','20261014 REVIEW','verification.json','PHYS-001/DICTIONARY.md','PHYS-003/PROOF.md'],'not_downloaded':['repo checkout','historical source packages','FIN-M0 trajectories'],'reconstructed_script':'run_campaign.py'})
    task_docs('RG-C01','model binding',current,'A7 constructed from frozen kernel and discrete Fourier rank-seven prescription. Centred norm error `%.3g`; rank %s; eigenvalues first seven %s. The kinetic rule is additionally *supplied*, not derived. The original source commit SHA is not verified.'%(np.max(np.abs(A.sum(axis=1))),np.linalg.matrix_rank(A),eig.tolist()),'A7 has no independent physical source. This is not a verified Git commit checkout, and no claims are imported from historical PASS labels.')
    static={};rows=[]
    for name,Am in models.items():
        for gain in [G,0.0]:
            F=marginal_F3(Am,gain);F2,R,decomp=three_decompose(F)
            p=prob(F);p2=prob(F2)
            oscR=osc(R);bound=math.tanh(oscR/4)
            s1,s2=sequential_direct_2(Am,gain)
            # exact four-to-three check against direct summation in entire 1728 states
            dmax=0.
            baseline=F[0,0,0]-logsumexp([labelled_log((0,0,0,d),4,gain,Am) for d in range(Q)])
            for a,b,c in itertools.product(range(Q),repeat=3):
                direct=logsumexp([labelled_log((a,b,c,d),4,gain,Am) for d in range(Q)])
                dmax=max(dmax,abs((F[a,b,c]-direct)-baseline))
            row={'model':name,'gain':gain,'contrast':alt_contrast(F),'pair_trunc_TV':tv(p,p2),'pair_osc':oscR,'TV_bound':bound,'direct_to_sequential_marginal_max_centered_diff':float(np.max(np.abs((s1-s1[0,0])-(s2-s2[0,0])))),'direct_labelled_marginal_max_abs_mod_additive_constant':dmax,'three_component_norm':decomp['three_body_inf_norm'],'decomposition':decomp}
            static[(name,gain)]=row;rows.append(row)
    (ROOT/'RG-C02').mkdir(exist_ok=True)
    with (ROOT/'RG-C02/GENERATED_OPERATORS.json').open('w') as f:json.dump(rows,f,indent=2)
    fin=static[('FIN_A7',G)];negative=[static[(k,0.0)] for k in models]
    statout={'execution_status':'EXECUTED','scientific_verdict':'NUMERICAL_FINITE_HIGHER_BODY_WITNESS_AND_CONTROLLED_TRUNCATION','model_id':current['model_id'],'input_hashes':{},'assumptions':BASE,'derived_objects':fin,'supplied_objects':['A7','gain'],'fitted_objects':[],'arithmetic_status':'binary64 finite exhaustive enumeration, interval source check separately','error_budget':{'finite_enumeration_pair_TV_bound':fin['TV_bound'],'actual_pair_TV':fin['pair_trunc_TV']},'discrepancies':[],'dependencies':['RG-C01'],'next_action':'test reduction dynamics'}
    task_docs('RG-C02','exact finite static reduction',statout,'The 4→3 elimination enumerates **all 1728 retained labelled states**, direct 4→3 sums of 12 and direct/sequential 4→2 sums of 144 removed states. FIN third contrast: **%.12g**. Pair-only approximation TV **%.9g**, residual oscillation **%.9g**, universal bound tanh(osc/4) **%.9g**. The decomposition is uniform-reference ANOVA/Hoeffding; positive third term is an effective model term, not a source. The negative-gain check covers all three models.\n\nProof of TV bound: If r=F-Ftilde has range width ω, shift its constant so |r|≤ω/2. The likelihood ratio of normalized measures has bounded log-range ω; the maximal TV across positive reweightings with that log-range is tanh(ω/4). A constant residual has zero error. Applying any fixed Markov observation kernel contracts TV.\n\n**Numerical status:** floating extrema are exhaustive over the finite alphabet, not machine-certified exact-kernel interval extrema; source-input enclosure is checked separately.'%(fin['contrast'],fin['pair_trunc_TV'],fin['pair_osc'],fin['TV_bound']),'Neither an emergent physical force, nor a spatial RG, a fixed point, or uniqueness of A7 follows. The TV envelope is mathematically valid conditional on true ω; current computed ω is numerical until fully outward-certified.')
    write_md('RG-C02/TRUNCATION_REPORT.md','\n'.join(['# Pair truncation by model','', '|model|gain|third contrast|actual TV|nominal TV bound|','|---|---:|---:|---:|---:|']+['|%s|%.9g|%.9g|%.9g|%.9g|'%(r['model'],r['gain'],r['contrast'],r['pair_trunc_TV'],r['TV_bound']) for r in rows]))
    (ROOT/'RG-C04').mkdir(exist_ok=True)
    with (ROOT/'RG-C04/COUNTERMODEL_COMPARISON.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['model','gain','contrast','pair_trunc_TV','pair_osc','TV_bound','three_component_norm']);w.writeheader();w.writerows([{k:r[k] for k in w.fieldnames} for r in rows])
    dy={name:dynamics(Am,G,full=(name=='FIN_A7')) for name,Am in models.items()}
    dfin=dy['FIN_A7'];dynout={'execution_status':'EXECUTED','scientific_verdict':'NUMERICAL_NONLUMPABILITY_AND_EXPLICIT_MEMORY_WITNESS','model_id':current['model_id'],'input_hashes':{},'assumptions':BASE,'derived_objects':dfin,'supplied_objects':['per-copy unit-rate leave-one-out kinetics'], 'fitted_objects':[],'arithmetic_status':'binary64 sparse CTMC and Krylov exponential','error_budget':{'simulation_residual':'independently assessed stationarity and detailed balance; no formal Krylov enclosure'},'discrepancies':[],'dependencies':['RG-C01'],'next_action':'countermodel comparison'}
    task_docs('RG-C03','projected process and memory',dynout,'For M=3, the sparse 1728-state CTMC is constructed directly from per-label heat-bath updates. Its coarse two-label observation has 144 values; the stationary conditional-averaged generator and the separately constructed marginal Gibbs sampler are both compared with actual projection. Rates for hidden states (0,0,0) and (0,0,1) differ, so strong lumpability **fails** at nonzero gain. The equilibrium distribution is shared by the two coarse generators, but finite-time transition distributions generally differ.\n\nThe scalar memory K(t)=⟨v,exp(t QLQ)v⟩_π, v=QLPf, f=1_{Y=(0,0)}, is evaluated for a reversible quadratic channel. K(0)=%.9g, derivative at 0 =%.9g. A nonequilibrium hidden initial condition is recorded separately; it is not silently assumed stationary.'%(dfin['memory_quadratic_channel']['M0_coupled'],dfin['memory_quadratic_channel']['M1_derivative']),'A stationary marginal does not uniquely determine observed dynamics. K(t) here is a single equilibrium quadratic observable, not a physical temporal memory or a universal coarse-time law. Dynamics are in generator time units, not seconds.')
    write_json('RG-C03/PROCESS_MAP.json',dy)
    write_json('RG-C03/LUMPABILITY_WITNESS.json',dfin['strong_lumpability_witness'])
    write_json('RG-C03/MEMORY_REPORT.json',dfin['memory_quadratic_channel'])
    comp=[]
    for name in models:
        s=static[(name,G)];d=dy[name]
        comp.append({'model':name,'generated_third_contrast':s['contrast'],'pair_TV':s['pair_trunc_TV'],'pair_TV_bound':s['TV_bound'],'nonlumpability_max_rate_difference':d['strong_lumpability_witness']['maximum_single_visible_refresh_rate_difference'],'equilibrium_memory_M0':d['memory_quadratic_channel']['M0_coupled'],'short_window_TV':d['finite_window'][0]['TV_exact_vs_equilibrium_generator']})
    specific={'execution_status':'EXECUTED','scientific_verdict':'MODEL_CLASS_MECHANISM_NOT_UNIQUE_TO_FIN','model_id':'FIN_A7 / FULL_POTTS_TRACE / FLAT_P7_TRACE','input_hashes':{},'assumptions':BASE,'derived_objects':comp,'supplied_objects':['three trace-matched matrices, identical gain and kinetics'],'fitted_objects':[],'arithmetic_status':'finite double precision','error_budget':{'source_matrix':'same model-class reconstruction'},'discrepancies':[],'dependencies':['RG-C02','RG-C03'],'next_action':'freeze computational holdout'}
    task_docs('RG-C04','trace-matched countermodels',specific,'The two controls are trace matched and use identical elimination, gain, observation and heat-bath kinetics. The generated higher-body contrast, strong-lumpability defect and/or memory channel are compared across all three. These effects are class-level unless absent from both controls. Full tabulation: COUNTERMODEL_COMPARISON.csv, RG-C03/PROCESS_MAP.json.','This is three finite comparisons, not a theorem of genericity or universal FIN specificity. A7 spectrum and its physical source remain free premises.')
    write_md('RG-C04/SPECIFICITY_VERDICT.md','# Scoped verdict\n\nThe mechanism of hidden-variable-induced memory and higher-body interactions is not automatically FIN-specific. Numerical contrasts vary with the trace-matched interaction matrix; cross-model differences can be reported but do not select nature’s A7. Compare the CSV and PROCESS_MAP before claiming a unique effect.\n')
    # HOLDOUT: freeze before running independent FOUR-LABEL outcome check
    test={'model':'FIN_A7','gain':G,'preparation':'stationary 4-copy Gibbs','observable':'1{first retained label == third retained label}','prediction_source':'3-copy 2-body ANOVA truncation only','held_out_reference':'independent direct 4-label sum, not the marginal_F3 routine','tolerance_rule':'absolute discrepancy <= nominal pair TV bound + 1e-10 numerical slack','controls':['FULL_POTTS_TRACE','FLAT_P7_TRACE'],'status':'frozen before heldout reference evaluations in this execution script','calibration':'no fitted parameters'}
    write_json('RG-C05/PREREGISTRATION.json',test)
    pred=[]
    for name,Am in models.items():
        F=marginal_F3(Am,G);F2,R,decomp=three_decompose(F)
        qpair=prob(F2); predicted=float(sum(qpair[i,j,i] for i in range(Q) for j in range(Q)))
        # direct reference only here, separate state enumeration 12^4, no reduction routines
        sums=0.;event=0.
        for state in itertools.product(range(Q),repeat=4):
            E=labelled_log(state,4,G,Am);weight=math.exp(E);sums+=weight
            if state[0]==state[2]: event+=weight
        measured=event/sums
        gap=abs(predicted-measured);bound=static[(name,G)]['TV_bound']
        pred.append({'model':name,'pair_trunc_prediction':predicted,'parent_direct_reference':measured,'abs_discrepancy':gap,'nominal_bound':bound,'passed_nominal_bound_plus_1e-10':gap<=bound+1e-10,'validation_kind':'synthetic finite-model computational holdout'})
    write_json('RG-C05/PREDICTIONS.json',{'frozen_contract':test,'predictions_and_reference':pred})
    validation={'execution_status':'EXECUTED','scientific_verdict':'COMPUTATIONAL_HOLDOUT_NOT_PHYSICAL_VALIDATION','model_id':current['model_id'],'input_hashes':{},'assumptions':BASE,'derived_objects':pred,'supplied_objects':['model parent enumeration'],'fitted_objects':[],'arithmetic_status':'binary64 enumeration, tolerance includes 1e-10 numerical slack','error_budget':{'truncation_bound_each':[r['nominal_bound'] for r in pred],'numeric_allowance':1e-10,'input_A7':'same rounded matrix in prediction/reference; exact-kernel external validation pending','memory':'not used, static holdout','readout':'exact event indicator'},'discrepancies':[],'dependencies':['RG-C02','RG-C04'],'next_action':'architect RG-C06 checkpoint (not decided by executor)'}
    task_docs('RG-C05','preregistered static computational holdout',validation,'The preregistered event was frozen in PREREGISTRATION.json before computing the independently enumerated 12^4 reference. Pair-truncation event predictions and exhaustive parent event expectations are compared with the same predeclared TV bound. This is a synthetic model-consistency holdout, not external empirical validation.','No claim of predictive success in the physical world, blind new-g test, or new independent microscopic source. A full exact-kernel interval validation would be required to call numerical bounds formally certified.')
    write_json('RG-C05/VALIDATION.json',validation)
    write_json('CLAIM_REGISTER.json',{'RG-C01':current['scientific_verdict'],'RG-C02':statout['scientific_verdict'],'RG-C03':dynout['scientific_verdict'],'RG-C04':specific['scientific_verdict'],'RG-C05':validation['scientific_verdict'],'RG-C06':'NOT_EXECUTED_ARCHITECT_DECISION','RG-M00':'OPEN_INPUT_NO_EXECUTABLE_MODEL','source_checkout':'NOT_AVAILABLE','Git_commit_SHA_latest':'NOT_VERIFIED'})
    write_md('RG-M00_OPEN_INPUT.md','# FIN-M0 not executed\n\nNo versioned simulator, controller, initial ensemble, trajectories or numerical specifications were supplied in the attached plan or made accessible in the code environment. Historical Kaggle/SOC/soliton claims are not replayed.\n')
    write_md('RG-C06_CHECKPOINT_FOR_ARCHITECT.md','# RG-C06 — architectural decision reserved\n\nRG-C01–RG-C05 export finite model evidence. The generator projection is non-Markov in general, so do not replace it with an exact two-label CTMC; retain hidden/memory information for dynamic predictions. The scientific architect must select a single branch in the separate RG-C06 decision. A source for A7 and any measurement bridge remain unprovided.\n')
    write_md('EVIDENCE_LEDGER.md','# Source/evidence ledger\n\n- The provided research plan is locally available and its SHA256 is recorded in root `SOURCE_ACCESS.md`.\n- Remote browsing: public `AGENTS.md`, methodology review, `verification.json`, PHYS-001 `DICTIONARY.md`, PHYS-002 `safe_seed.py`, PHYS-003 `PROOF.md` were inspected through GitHub Raw. They could not be downloaded into the execution runtime for byte-for-byte commit-specific hashing.\n- Code in this package independently implements the published safe-seed formula; no claim of source identity is made.\n- All displayed numerical data were generated in this package, not copied from an unreplayed PASS report.\n- Historical PHYS-013 calibration and FIN-M0 results are not considered reopened tests.\n')
    write_md('HANDOFF.md','# FIN generative campaign — execution handoff (2026-10-14)\n\n## Scope / status\n\nRG-C01 to RG-C05 are executed against an independently reconstructed q=12 rank-seven finite model using published formulas, not a cloned/verified checkout. RG-C06 remains reserved for the architect, RG-M00 is blocked on sources.\n\n## Results\n\n- RG-C01: supplied A7 Fourier construction, gain and per-copy leave-one-out CTMC frozen; source checkout SHA unavailable.\n- RG-C02: FIN third contrast %.12g (earlier review ≈0.27127058458855), actual pair TV %.8g; mathematically valid tanh(osc/4) nominal certificate %.8g. Source interval checked separately where applicable.\n- RG-C03: nonzero strong-lumpability defect %.8g, K(0)=%.8g and explicit finite-time deviations of observed dynamics from averaged Markov model.\n- RG-C04: trace-matched Potts and flat rank-seven numerical controls; see CSV.\n- RG-C05: preregistered static computational holdout: FIN |error| %.8g vs pair TV envelope %.8g.\n\n## Blockers and next action\n\nThe latest commit SHA and full original source/input hashes could not be confirmed: raw files were viewable, `git clone` and GitHub commit history were unavailable. Independently sourced physical A7, FIN-M0 simulator/trajectories, an exact-kernel outward interval certificate for the full F3 oscillation and external experimental data remain absent. Do not promote any of these. RG-C06 requires an architect to decide whether the next branch retains memory, restores FIN-M0, or targets an engineered readout.\n\n## Replay\n\n`python run_campaign.py` recreates numerical outputs; `python independent_verify.py` independently checks claimed invariants. Versions and checksums in `ENVIRONMENT.json`, `INPUTS.sha256`, and `MANIFEST.sha256`.\n'%(fin['contrast'],fin['pair_trunc_TV'],fin['TV_bound'],dfin['strong_lumpability_witness']['maximum_single_visible_refresh_rate_difference'],dfin['memory_quadratic_channel']['M0_coupled'],pred[0]['abs_discrepancy'],pred[0]['nominal_bound']))
    import scipy,mpmath
    write_json('ENVIRONMENT.json',{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mpmath.__version__,'runtime_seconds':time.time()-started})
    print('FINISHED',time.time()-started,'C02',fin['contrast'],fin['pair_trunc_TV'],fin['TV_bound'],'C03 memory',dfin['memory_quadratic_channel']['M0_coupled'],'holdout',pred[0])
if __name__=='__main__':main()
