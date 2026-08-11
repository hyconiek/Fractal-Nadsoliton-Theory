#!/usr/bin/env python3
"""FIN ST118--ST129: operational equivalence, certified design, and typed tests."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from mpmath.libmp.libmpf import ComplexResult
from scipy.optimize import root

from fin_st01_st15_research import N, strict_operator
from fin_st28_st45_research import dyadic_lift
from fin_st106_st117_research import plus_minus
from fin_st120_minimal_replay import run_replay


ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"FIN_ST118_ST129_Results.json"
SUMMARY=ROOT/"FIN_ST118_ST129_Summary.csv"
PACK119=ROOT/"FIN_ST119_Circulant_Lift_Polyhedron.json"
CERT124=ROOT/"FIN_ST124_Continuous_Beta_Interval_Certificate.json"
PACK129=ROOT/"FIN_ST129_Minimax_Fine_Observable.json"
FIG_DIR=ROOT/"FIN_ST118_ST129_Figures"
SEED=20260821


def native(x:Any)->Any:
    if isinstance(x,dict): return {str(k):native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [native(v) for v in x]
    if isinstance(x,np.ndarray): return native(x.tolist())
    if isinstance(x,(np.floating,np.integer)): return x.item()
    return x


def sha(path:Path)->str: return hashlib.sha256(path.read_bytes()).hexdigest()


def ib(x):
    if isinstance(x,tuple): return mp.iv.mpf([str(x[0]),str(x[1])])
    if isinstance(x,(float,np.floating)): return mp.iv.mpf([str(float(np.nextafter(x,-np.inf))),str(float(np.nextafter(x,np.inf)))])
    return mp.iv.mpf(str(x))


def bounds(x): return float(np.nextafter(float(x.a),-np.inf)),float(np.nextafter(float(x.b),np.inf))


def interval_product(a0,a1,b0,b1):
    vals=[a0*b0,a0*b1,a1*b0,a1*b1]
    return float(np.nextafter(min(vals),-np.inf)),float(np.nextafter(max(vals),np.inf))


def interval_matvec(alo,ahi,xlo,xhi):
    lo=np.zeros(alo.shape[0]); hi=np.zeros(alo.shape[0])
    for i in range(alo.shape[0]):
        l=h=0.0
        for j in range(alo.shape[1]):
            q0,q1=interval_product(alo[i,j],ahi[i,j],xlo[j],xhi[j])
            l=np.nextafter(l+q0,-np.inf); h=np.nextafter(h+q1,np.inf)
        lo[i],hi[i]=l,h
    return lo,hi


def interval_left(r,alo,ahi):
    lo=np.zeros((r.shape[0],alo.shape[1])); hi=np.zeros_like(lo)
    for i in range(r.shape[0]):
        for j in range(alo.shape[1]):
            l=h=0.0
            for k in range(r.shape[1]):
                q0,q1=interval_product(r[i,k],r[i,k],alo[k,j],ahi[k,j])
                l=np.nextafter(l+q0,-np.inf); h=np.nextafter(h+q1,np.inf)
            lo[i,j],hi[i,j]=l,h
    return lo,hi


def st118_operational_equivalence()->dict:
    return {
        "program":"ST118","object":"Universal Fine-Blind Operational Equivalence",
        "group":"G_f={I_c direct-sum U_f: U_f in U(12)}",
        "fixed_algebra":"M12 direct-sum C","fixed_algebra_complex_dimension":145,
        "haar_expectation":"E(X)=Pc X Pc + Tr(Pf X Pf) Pf/12",
        "theorem":(
            "If two observables are operationally equivalent exactly when all states and instruments invariant under arbitrary "
            "fine-sector basis changes give the same records, the separating quotient is the fixed-point algebra of G_f. "
            "Schur's lemma makes this algebra M12 direct-sum C, and Haar twirling gives the unique trace-preserving expectation."
        ),
        "necessity_test":(
            "Removing arbitrary fine U(12) invariance and retaining only layer swap restores M12 direct-sum M12; removing "
            "fine blindness entirely restores M24.  Thus the equivalence axiom is sufficient and minimal relative to those two weakenings."
        ),
        "status":"proven_universal_property_conditional_on_explicit_fine_gauge_equivalence",
        "boundary":"The operational equivalence is an added gauge/blindness axiom, not selected by strict A or a laboratory apparatus.",
    }


def circulant_eigenvalues(first:np.ndarray)->np.ndarray:
    return np.real(np.fft.fft(first))


def st119_circulant_polyhedron(a:np.ndarray)->dict:
    weights=np.array([-a[0,d] for d in range(1,7)]); s=float(a[0,0])
    vertices=[]
    for signs in itertools.product([-1,1],repeat=6):
        b=np.r_[s,weights*np.array(signs)]
        first=np.r_[b[:7],b[5:0:-1]]
        vertices.append({"signs":list(signs),"parameters":b.tolist(),"minimum_eigenvalue":float(np.min(circulant_eigenvalues(first)))})
    packet={
        "radial_parameters":"b0,b1,...,b6","vertex_count":len(vertices),"recession_extreme_ray":"positive b0 direction",
        "vertices":vertices,"minimum_vertex_eigenvalue":min(v["minimum_eigenvalue"] for v in vertices),
        "theorem":(
            "In the real symmetric circulant section, graph bounds give |b_d|<=w_d for d=1..6 and b0>=s. "
            "Because s=2 sum_(d=1)^5 w_d+w_6, every such matrix is symmetric diagonally dominant with nonnegative diagonal "
            "and hence PSD.  The section is exactly [s,infinity) times the six-dimensional box product[-w_d,w_d]. "
            "It has 2^6=64 vertices at b0=s and b_d=plus/minus w_d, plus one recession ray in b0."
        ),
        "scope":"Exact structural classification conditional on the declared positive strict weights; displayed eigenvalues are binary64 checks.",
    }
    PACK119.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":"ST119","object":"Complete Circulant Section of the Lift Spectrahedron","packet_file":PACK119.name,
            "packet_sha256":sha(PACK119),**packet,"status":"proven_exact_circulant_polyhedral_classification",
            "boundary":"This seven-parameter section does not enumerate all noncirculant faces of the 78-dimensional fiber."}


def st120_replay()->dict:
    replay=run_replay(True)
    path=ROOT/"FIN_ST120_Minimal_Transcendental_Replay_Certificate.json"
    return {"program":"ST120","object":"Dependency-Minimal Independent ST108 Replay","certificate_file":path.name,
            "certificate_sha256":sha(path),**replay,
            "status":"proven_independent_two_dependency_interval_replay" if replay["accepted"] else "blocked_replay_failed",
            "boundary":"The numerical root center remains an input; this is not a standard-library-only or proof-assistant replay."}


def st121_code_normal_form()->dict:
    return {
        "program":"ST121","object":"Normal Form of Every Maximal Pinching-Recoverable Code",
        "classification":{
            "alpha_0":"the full fine sector, up to a fine unitary",
            "alpha_1":"the full coarse sector, up to a coarse unitary",
            "0_lt_alpha_lt_1":"graph code V psi=sqrt(alpha)P Uc psi+sqrt(1-alpha)F Uf psi",
            "equivalence_invariant":"alpha, after independent sector unitaries",
        },
        "theorem":(
            "For a twelve-dimensional correctable code Q, Knill--Laflamme gives QPcQ=alpha Q.  When 0<alpha<1, "
            "the normalized coarse and fine restrictions are unitary maps onto their rank-twelve sectors, yielding the graph-code normal form. "
            "At alpha=0 or 1 the code is one whole sector.  Conversely every displayed code is correctable."
        ),
        "status":"proven_complete_maximal_code_normal_form",
        "boundary":"The continuous invariant alpha and sector unitaries are not selected by strict FIN.",
    }


def st122_invariant_state_noselector(a:np.ndarray)->dict:
    rows=[]
    for beta in [0,.2,.8,2,8]:
        vals,vecs=np.linalg.eigh(a); rho=vecs@np.diag(np.exp(-beta*vals))@vecs.T; rho/=np.trace(rho)
        diag=np.diag(rho)
        rows.append({"beta":beta,"maximum_branch_probability_deviation":float(np.max(abs(diag-1/N))),
                     "reflection_commutator":float(np.linalg.norm(rho-np.flip(np.flip(rho,0),1)))})
    return {
        "program":"ST122","object":"Strict Functional-Calculus State Nonselection Theorem","rows":rows,
        "theorem":(
            "Every state rho=f(A)/Tr f(A) obtained from strict functional calculus is invariant under the transitive C12 action "
            "and reflection.  Any covariant twelve-branch vertex POVM therefore has uniform outcome probabilities.  More generally, "
            "an invariant state paired with a covariant transitive instrument cannot generate a nonuniform branch law."
        ),
        "status":"proven_no_nonuniform_selector_in_strict_invariant_state_class",
        "boundary":"Noninvariant initial/boundary states or instruments evade the theorem only as newly supplied selector sources; QW-2191 remains open.",
    }


def st123_minimal_active_extension()->dict:
    gamma,kappa=.35,1.0
    rows=[]
    for g in np.linspace(0,.9,19):
        growth=float(g-gamma); growth=0.0 if abs(growth)<1e-14 else growth
        amp=math.sqrt(max(growth,0)/kappa)
        rows.append({"pump_g":float(g),"linear_growth_rate":growth,"stable_nonzero_amplitude":amp,
                     "pump_injection_at_attractor":float(g*amp*amp),"dissipation_at_attractor":float(gamma*amp*amp),
                     "saturation_at_attractor":float(kappa*amp**4)})
    return {
        "program":"ST123","object":"Minimal Scalar Pump and Saturation Extension","parameters":{"gamma":gamma,"kappa":kappa},"rows":rows,
        "equation":"x_dot=(g-gamma)x-kappa x^3+u, y=x",
        "theorem":(
            "For g<=gamma the zero state is stable and no self-amplified branch exists.  At g=gamma a pitchfork occurs; for g>gamma "
            "the zero state is unstable and the two equilibria plus/minus sqrt((g-gamma)/kappa) are stable.  The storage balance explicitly "
            "contains pump input g x^2, dissipation gamma x^2 and saturation kappa x^4.  Active gain therefore requires a sourced pump."
        ),
        "status":"proven_constructed_minimal_active_bifurcation_and_cost_balance",
        "boundary":"g, gamma, kappa, the clock and energetic units are added; this is not a strict FIN-derived information amplifier.",
    }


def spectrum_intervals(vals): return [ib(float(v)) for v in vals]


def partition_moments(beta,eigs):
    terms=[mp.iv.exp(-beta*x) for x in eigs]; z=sum(terms,ib(0))
    m1=sum((x*t for x,t in zip(eigs,terms)),ib(0))/z
    m2=sum((x*x*t for x,t in zip(eigs,terms)),ib(0))/z
    return z,m1,m2-m1*m1


def detected_probability_jet(beta,q,eigc,eigf,delta):
    zc,mc,vc=partition_moments(beta,eigc); zf,mf,vf=partition_moments(beta,eigf)
    lr=mp.iv.exp(-2*beta*ib(q))*zf/zc; p=1/(1+lr)
    lp=-2*ib(q)-mf+mc; lpp=vf-vc
    dp=-p*(1-p)*lp; ddp=p*(1-p)*((1-2*p)*lp*lp-lpp)
    c=1-2*ib(delta)
    return ib(delta)+c*p,c*dp,c*ddp


def chernoff_stationary(beta,s,eigc,eigf,delta=.05):
    a,ap,app=detected_probability_jet(beta,.2,eigc,eigf,delta)
    b,bp,bpp=detected_probability_jet(beta,.7,eigc,eigf,delta)
    Bs=ib(0); Bb=ib(0); Bss=ib(0); Bsb=ib(0); Bbb=ib(0); B=ib(0)
    for x,xp,xpp,y,yp,ypp in [(a,ap,app,b,bp,bpp),(1-a,-ap,-app,1-b,-bp,-bpp)]:
        lx=mp.iv.log(x); ly=mp.iv.log(y); term=mp.iv.exp(s*lx+(1-s)*ly)
        hs=lx-ly; hb=s*xp/x+(1-s)*yp/y
        hsb=xp/x-yp/y
        hbb=s*(xpp/x-(xp/x)**2)+(1-s)*(ypp/y-(yp/y)**2)
        B+=term; Bs+=term*hs; Bb+=term*hb; Bss+=term*hs**2; Bsb+=term*(hb*hs+hsb); Bbb+=term*(hb**2+hbb)
    return [Bs,Bb],[[Bss,Bsb],[Bsb,Bbb]],B


def point_system(x,eigc_float,eigf_float):
    beta,s=x
    def jet(q):
        tc=np.exp(-beta*eigc_float); tf=np.exp(-beta*eigf_float); zc=tc.sum(); zf=tf.sum()
        mc=(eigc_float*tc).sum()/zc; mf=(eigf_float*tf).sum()/zf
        vc=(eigc_float**2*tc).sum()/zc-mc**2; vf=(eigf_float**2*tf).sum()/zf-mf**2
        lr=np.exp(-2*beta*q)*zf/zc; p=1/(1+lr); lp=-2*q-mf+mc; lpp=vf-vc
        dp=-p*(1-p)*lp; ddp=p*(1-p)*((1-2*p)*lp**2-lpp); c=.9
        return .05+c*p,c*dp,c*ddp
    a,ap,app=jet(.2); b,bp,bpp=jet(.7); Bs=Bb=0.0
    for xx,xp,xpp,yy,yp,ypp in [(a,ap,app,b,bp,bpp),(1-a,-ap,-app,1-b,-bp,-bpp)]:
        term=xx**s*yy**(1-s); Bs+=term*math.log(xx/yy); Bb+=term*(s*xp/xx+(1-s)*yp/yy)
    return np.array([Bs,Bb])


def sstar_interval(beta,eigc,eigf):
    a,ap,_=detected_probability_jet(beta,.2,eigc,eigf,.05); b,bp,_=detected_probability_jet(beta,.7,eigc,eigf,.05)
    lp=mp.iv.log(a/b); ln=mp.iv.log((1-a)/(1-b)); ratio=ln/(-lp)
    s=(mp.iv.log(ratio)-mp.iv.log(b/(1-b)))/(mp.iv.log(a/(1-a))-mp.iv.log(b/(1-b)))
    return s,(a,ap,b,bp)


def envelope_bbeta(beta,eigc,eigf):
    s,(a,ap,b,bp)=sstar_interval(beta,eigc,eigf); total=ib(0)
    for x,xp,y,yp in [(a,ap,b,bp),(1-a,-ap,1-b,-bp)]:
        term=x**s*y**(1-s); total += term*(s*xp/x+(1-s)*yp/y)
    return total,s


def st124_interval_beta(a:np.ndarray)->dict:
    p,f=plus_minus(); b0=f.T@dyadic_lift(a,0)@f
    ec=np.linalg.eigvalsh(a); ef=np.linalg.eigvalsh(b0); eigc=spectrum_intervals(ec); eigf=spectrum_intervals(ef)
    center=root(lambda x:point_system(x,ec,ef),np.array([2.1828,.5398]),tol=1e-12).x
    mp.iv.dps=60; radii=[1e-6,3e-7,1e-7,3e-8]; attempts=[]; accepted=None
    # Krawczyk for (beta,s), system (B_s,B_beta).
    f0,j0,_=chernoff_stationary(ib(float(center[0])),ib(float(center[1])),eigc,eigf)
    flo=np.array([bounds(x)[0] for x in f0]); fhi=np.array([bounds(x)[1] for x in f0])
    jlo=np.array([[bounds(x)[0] for x in row] for row in j0]); jhi=np.array([[bounds(x)[1] for x in row] for row in j0]); c=np.linalg.inv((jlo+jhi)/2)
    for rad in radii:
        fi,ji,_=chernoff_stationary(ib((center[0]-rad,center[0]+rad)),ib((center[1]-rad,center[1]+rad)),eigc,eigf)
        jl=np.array([[bounds(x)[0] for x in row] for row in ji]); jh=np.array([[bounds(x)[1] for x in row] for row in ji])
        cglo,cghi=interval_matvec(c,c,flo,fhi); ylo=np.nextafter(center-cghi,-np.inf); yhi=np.nextafter(center-cglo,np.inf)
        cjlo,cjhi=interval_left(c,jl,jh); mlo,mhi=-cjhi,-cjlo
        for i in range(2): mlo[i,i]=np.nextafter(mlo[i,i]+1,-np.inf); mhi[i,i]=np.nextafter(mhi[i,i]+1,np.inf)
        mdlo,mdhi=interval_matvec(mlo,mhi,np.full(2,-rad),np.full(2,rad)); klo=ylo+mdlo; khi=yhi+mdhi
        margin=float(min(np.min(klo-(center-rad)),np.min((center+rad)-khi)))
        row={"radius":rad,"margin":margin,"included":margin>0}; attempts.append(row)
        if accepted is None and margin>0: accepted=row
    # Global derivative-sign cover outside the accepted local beta interval.
    root_rad=accepted["radius"] if accepted else 1e-6; left_end=center[0]-root_rad; right_start=center[0]+root_rad
    cover=[]; unresolved=[]
    def certify(lo,hi,side,depth=0):
        try:
            val,_=envelope_bbeta(ib((lo,hi)),eigc,eigf); vl,vh=bounds(val)
        except (ValueError, ZeroDivisionError, ComplexResult):
            if depth>=2 or hi-lo<2e-6:
                unresolved.append({"lo":lo,"hi":hi,"B_beta_interval":None,"side":side,"reason":"ill_conditioned_s_formula"}); return
            mid=(lo+hi)/2; certify(lo,mid,side,depth+1); certify(mid,hi,side,depth+1); return
        ok=(vh<0) if side=="left" else (vl>0)
        if ok: cover.append({"lo":lo,"hi":hi,"B_beta_interval":[vl,vh],"side":side}); return
        if depth>=2 or hi-lo<2e-6: unresolved.append({"lo":lo,"hi":hi,"B_beta_interval":[vl,vh],"side":side}); return
        mid=(lo+hi)/2; certify(lo,mid,side,depth+1); certify(mid,hi,side,depth+1)
    grid=np.linspace(.05,left_end,61)
    for lo,hi in zip(grid[:-1],grid[1:]): certify(float(lo),float(hi),"left")
    grid=np.linspace(right_start,5.0,81)
    for lo,hi in zip(grid[:-1],grid[1:]): certify(float(lo),float(hi),"right")
    # Independent dense floating audit; evidence only, never used as the interval proof gate.
    def point_s_and_bbeta(beta):
        def prob(q):
            tc=np.exp(-beta*ec); tf=np.exp(-beta*ef); zc=tc.sum(); zf=tf.sum()
            mc=(ec*tc).sum()/zc; mf=(ef*tf).sum()/zf
            lr=np.exp(-2*beta*q)*zf/zc; p=1/(1+lr); dp=-p*(1-p)*(-2*q-mf+mc)
            return .05+.9*p,.9*dp
        aa,aap=prob(.2); bb,bbp=prob(.7); lp=math.log(aa/bb); ln=math.log((1-aa)/(1-bb))
        ss=(math.log(ln/(-lp))-math.log(bb/(1-bb)))/(math.log(aa/(1-aa))-math.log(bb/(1-bb)))
        return ss,float(point_system(np.array([beta,ss]),ec,ef)[1])
    sampled=[]
    for beta in np.linspace(.05,5,2001):
        ss,bbeta=point_s_and_bbeta(float(beta)); sampled.append((float(beta),ss,bbeta))
    sign_changes=sum(sampled[i][2]*sampled[i+1][2]<0 for i in range(len(sampled)-1))
    _,_,bpoint=chernoff_stationary(ib(float(center[0])),ib(float(center[1])),eigc,eigf)
    cert={"domain":[.05,5.0],"stationary_center":{"beta":float(center[0]),"s":float(center[1])},"krawczyk_attempts":attempts,"accepted":accepted,
          "derivative_sign_boxes":len(cover),"unresolved_boxes":unresolved,"Chernoff_information_interval_at_center":list(bounds(-mp.iv.log(bpoint))),
          "sampled_beta_points":len(sampled),"sampled_B_beta_sign_changes":sign_changes,
          "proof_method":"2D interval Krawczyk for B_s=B_beta=0 plus global interval sign cover of envelope B_beta outside the root box",
          "scope":"Frozen binary64 strict/dyadic spectra, q0=0.2, q1=0.7, symmetric detector flip delta=0.05, beta in [0.05,5]."}
    CERT124.write_text(json.dumps(native(cert),indent=2,sort_keys=True),encoding="utf-8")
    passed=accepted is not None and not unresolved
    return {"program":"ST124","object":"Continuous-Beta Chernoff Optimum Certificate","certificate_file":CERT124.name,"certificate_sha256":sha(CERT124),**cert,
            "status":"proven_unique_global_beta_optimum_in_declared_binary64_model" if passed else "strong_local_certificate_global_cover_incomplete",
            "boundary":"No physical inverse temperature, energy scale, preparation, or detector calibration is derived."}


def st125_mixed_refinement()->dict:
    degrees=[2,3,5,4,3]
    rows=[]; h=-1
    for n,m in enumerate(degrees):
        new=h**m; rows.append({"level":n,"degree":m,"input_h":h,"untwisted_output_h":new,"sign_survives":new==h}); h=new
    return {"program":"ST125","object":"Mixed-Degree Z2 Refinement Classification","example_rows":rows,
            "theorem":(
                "For a degree-m refinement, Z2 holonomy pulls back as h maps to h^m.  Odd degrees preserve h; any even degree sends both "
                "classes to +1, after which untwisted refinements remain periodic.  In a mixed tower, an initial antiperiodic class survives "
                "exactly until the first even degree.  Twisted recursion h_(n+1)=t_n h_n^(m_n) again makes every recovery after an even step new data."
            ),"status":"proven_mixed_refinement_parity_theorem","boundary":"The degree sequence and twists are supplied; no strict spin source follows."}


def st126_gauge_net_locality()->dict:
    return {"program":"ST126","object":"Gauge-Net Alternative to Vertex Locality",
            "net":"A(O)=P M_O P* direct-sum C Pf for coarse regions O; fine U(12) is a global gauge action",
            "theorem":(
                "The gauge-fixed local net is isotone, and disjoint coarse matrix algebras commute when the chosen coarse regions do. "
                "Its global algebra is M12 direct-sum C.  However, removing the fine gauge action while retaining the same coarse net gives "
                "M12 direct-sum M12, and adjoining vertex-resolving fine projectors generates M24.  Locality alone is insufficient; gauge equivalence does all selecting work."
            ),"status":"proven_conditional_gauge_net_construction_and_locality_no_go",
            "boundary":"The construction reorganizes the ST118 axiom and does not derive it from strict dynamics or causality."}


def st127_correlated_nuisance()->dict:
    p0,p1=.72,.94
    one_event_error=.5*(1-abs(p1-p0))
    rows=[]
    for n in [1,2,5,10,50,100,1000]: rows.append({"events":n,"iid_error_upper_heuristic":float(.5*math.exp(-n*.02)),"fully_correlated_exact_error":one_event_error})
    return {"program":"ST127","object":"Correlated-Event Error-Exponent Impossibility",
            "example_marginals":{"p0":p0,"p1":p1},"rows":rows,
            "theorem":(
                "One-event marginals and detector-error bounds do not imply a positive asymptotic discrimination exponent. "
                "Let one Bernoulli variable Y be drawn under each hypothesis and report X_1=...=X_n=Y.  Every event has the declared marginal, "
                "but the n-event likelihood contains only one independent bit, so the optimal Bayes error equals the one-event value for every n."
            ),"status":"proven_no_positive_count_exponent_under_arbitrary_temporal_correlation",
            "boundary":"ST115 remains valid per joint event; finite-count scaling additionally requires independence or a certified mixing bound."}


def st128_typed_thermal_class(a:np.ndarray)->dict:
    beta=.8; vals=np.linalg.eigvalsh(a); unique=[]
    for x in vals:
        if not unique or abs(x-unique[-1][0])>1e-9: unique.append([float(x),1])
        else: unique[-1][1]+=1
    gamma=np.exp(-beta*vals); gamma/=gamma.sum()
    return {"program":"ST128","object":"Exact Separation for a Declared Zero-Hamiltonian Bath Class",
            "declared_class":"finite baths of arbitrary finite dimension with H_B=0 and exact [U,H_S tensor I_B]=0",
            "system_energy_multiplicities":[u[1] for u in unique],"beta":beta,"gibbs_min_probability":float(gamma.min()),
            "witness_channel":"Phi(rho)=gamma Tr(rho)",
            "theorem":(
                "The reset-to-Gibbs channel is CPTP, Gibbs preserving, and time-translation covariant.  Any dilation in the declared H_B=0 class "
                "commutes with every system spectral projector and therefore preserves each system energy-block population for every input. "
                "The reset channel changes those populations unless the input already has Gibbs populations, so it lies outside this declared thermal subclass."
            ),"status":"proven_CGP_outside_declared_zero_energy_bath_subclass",
            "boundary":"This does not separate CGP from standard thermal operations with energy-bearing baths, catalysts, approximations, or unbounded ancillary classes."}


def st129_minimax_detector()->dict:
    rng=np.random.default_rng(SEED+129); trials=[]
    for _ in range(1000):
        h=rng.dirichlet(np.ones(N)); trials.append(float(np.min(h)))
    packet={"adversary":"t_i>=0, sum_i t_i=1 in Delta L=F diag(t) F*","detector_budget":"H=F diag(h)F*, h_i>=0, sum_i h_i=1",
            "optimal_h":[1/N]*N,"minimax_value":1/N,"best_random_trial":max(trials),
            "theorem":(
                "The detector response is h dot t, whose minimum over the adversarial simplex is min_i h_i.  Under sum h_i=1, "
                "min_i h_i<=1/12, with equality only for uniform h.  Thus H=Pf/12 is the unique diagonal minimax detector and guarantees response 1/12. "
                "Every coarse observable has zero response; every single-site detector has worst-case response zero."
            )}
    PACK129.write_text(json.dumps(native(packet),indent=2,sort_keys=True),encoding="utf-8")
    return {"program":"ST129","object":"Minimax Fine Observable against the Adversarial Lift Simplex","packet_file":PACK129.name,"packet_sha256":sha(PACK129),**packet,
            "status":"proven_exact_minimax_diagonal_fine_detector","boundary":"The normalized adversarial cone and detector budget are supplied; no laboratory POVM or noise model is derived."}


def make_figures(d):
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4)); ax.bar(["fine-blind","swap","vertex"],[145,288,576]); ax.set(ylabel="complex dimension",title="ST118/ST126: equivalence, not locality, selects the algebra"); fig.tight_layout(); fig.savefig(FIG_DIR/"st118_equivalence_algebras.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); vals=[v["minimum_eigenvalue"] for v in d["ST119"]["vertices"]]; ax.hist(vals,bins=15); ax.set(xlabel="minimum vertex eigenvalue",ylabel="count",title="ST119: 64 exact circulant vertices"); fig.tight_layout(); fig.savefig(FIG_DIR/"st119_vertices.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=d["ST120"]["attempts"]; ax.semilogx([r["radius"] for r in rows],[r["minimum_strict_inclusion_margin"] for r in rows],"o-"); ax.axhline(0,color="black"); ax.set(xlabel="box radius",ylabel="inclusion margin",title="ST120 independent replay"); fig.tight_layout(); fig.savefig(FIG_DIR/"st120_replay.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=d["ST123"]["rows"]; ax.plot([r["pump_g"] for r in rows],[r["stable_nonzero_amplitude"] for r in rows]); ax.axvline(d["ST123"]["parameters"]["gamma"],color="red",ls="--"); ax.set(xlabel="pump g",ylabel="stable amplitude",title="ST123 pump threshold and saturated branches"); fig.tight_layout(); fig.savefig(FIG_DIR/"st123_active_threshold.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); c=d["ST124"]; ax.bar(["beta","s"],[c["stationary_center"]["beta"],c["stationary_center"]["s"]]); ax.set(title="ST124 certified stationary design",ylabel="value"); fig.tight_layout(); fig.savefig(FIG_DIR/"st124_beta_certificate.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=d["ST125"]["example_rows"]; ax.step(range(len(rows)),[r["untwisted_output_h"] for r in rows],where="mid"); ax.set(ylim=(-1.2,1.2),xlabel="refinement step",ylabel="holonomy",title="ST125 first even degree erases antiperiodicity"); fig.tight_layout(); fig.savefig(FIG_DIR/"st125_mixed_refinement.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); rows=d["ST127"]["rows"]; ax.semilogx([r["events"] for r in rows],[r["fully_correlated_exact_error"] for r in rows],"o-",label="fully correlated"); ax.semilogx([r["events"] for r in rows],[r["iid_error_upper_heuristic"] for r in rows],"--",label="i.i.d. heuristic"); ax.set(xlabel="events",ylabel="error",title="ST127 correlations destroy count exponent"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR/"st127_correlated_nuisance.png",dpi=190); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4)); ax.bar(["single site","best random","uniform minimax"],[0,d["ST129"]["best_random_trial"],d["ST129"]["minimax_value"]]); ax.set(ylabel="guaranteed response",title="ST129 minimax fine detector"); fig.tight_layout(); fig.savefig(FIG_DIR/"st129_minimax_detector.png",dpi=190); plt.close(fig)


def main():
    _,a,_=strict_operator(); out={"metadata":{"programs":"ST118-ST129","date":"2026-08-11","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    out["ST118"]=st118_operational_equivalence(); out["ST119"]=st119_circulant_polyhedron(a); out["ST120"]=st120_replay(); out["ST121"]=st121_code_normal_form()
    out["ST122"]=st122_invariant_state_noselector(a); out["ST123"]=st123_minimal_active_extension(); out["ST124"]=st124_interval_beta(a); out["ST125"]=st125_mixed_refinement()
    out["ST126"]=st126_gauge_net_locality(); out["ST127"]=st127_correlated_nuisance(); out["ST128"]=st128_typed_thermal_class(a); out["ST129"]=st129_minimax_detector()
    out["recommended_next_programs"]=[
        {"id":"ST130","priority":1,"study":"seek a strict operational source for the fine U(12) gauge equivalence, or prove none exists in functional calculus"},
        {"id":"ST131","priority":2,"study":"enumerate a noncirculant low-rank face family of the full 78-dimensional lift spectrahedron"},
        {"id":"ST132","priority":3,"study":"replace the imported ST108 center in the minimal replay by an interval Newton center-isolation stage"},
        {"id":"ST133","priority":4,"study":"derive optimal recovery when coarse/fine pinching is noisy rather than exact"},
        {"id":"ST134","priority":5,"study":"classify symmetry-breaking states outside strict functional calculus by minimal tensor rank and source cost"},
        {"id":"ST135","priority":6,"study":"couple the scalar pump to strict modes and prove or obstruct orbital stability without claiming physical gain"},
        {"id":"ST136","priority":7,"study":"extend the beta certificate over interval detector flip and interval q hypotheses"},
        {"id":"ST137","priority":8,"study":"classify coherent twists for arbitrary finite refinement diagrams and path independence"},
        {"id":"ST138","priority":9,"study":"test whether a causal conditional-independence axiom selects the gauge net without naming the fine sector"},
        {"id":"ST139","priority":10,"study":"derive finite-count bounds under certified alpha-mixing or Markov detector correlations"},
        {"id":"ST140","priority":11,"study":"increase the typed bath from H_B=0 to one matched-gap qubit and solve exact dilation feasibility"},
        {"id":"ST141","priority":12,"study":"minimax-observable design over bounded non-diagonal lift faces and measurement noise"},
    ]
    out["central_verdict"]=("An explicit fine-sector gauge equivalence uniquely yields the desired observable algebra, but the equivalence remains an axiom and locality does not replace it. "
        "The circulant hidden-lift section is exactly a six-cube times one ray with 64 vertices.  The transcendental fold independently replays, maximal recoverable codes have a complete normal form, "
        "strict invariant states remain nonselecting, active gain requires a pump, and arbitrary temporal correlation destroys finite-count exponents.  A typed zero-energy bath subclass permits an exact CGP separation only in that restricted class.")
    out["epistemic_boundary"]="No strict fine-gauge source, QW-2191 closure, dimensional unit, apparatus, laboratory evidence, legacy-to-strict completion or role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h); w.writerow(["program","object","status"])
        for k in range(118,130): w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out); print(json.dumps({"results":RESULTS.name,"programs":12,"figures":8},indent=2))


if __name__=="__main__": main()
