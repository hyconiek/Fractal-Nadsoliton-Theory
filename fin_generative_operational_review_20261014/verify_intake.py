#!/usr/bin/env python3
"""Immutable RG-E intake: versioned inputs, exact N2 quotient, interval cover.

No imports from the new producer. Interval certificates assume mpmath.iv.
Dynamic matrix operations remain numerical. All new outputs are separate.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib.util
import itertools
import json
import math
import re
import zipfile
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.special import softmax

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "fin_generative_operational_campaign_20261014"
Q = 12


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def endpoints(z):
    out = []
    for s,m,e,_ in z._mpi_:
        n = -m if s else m
        out.append(Fraction(n << e,1) if e >= 0 else Fraction(n,1 << -e))
    return out


def decimal(x, upper=True, digits=15):
    n,d = x.numerator*10**digits,x.denominator
    k = -((-n)//d) if upper else n//d
    sign = "-" if k < 0 else ""
    a,b = divmod(abs(k),10**digits)
    return f"{sign}{a}.{b:0{digits}d}"


def box(z):
    a,b = endpoints(z)
    return [decimal(a,False),decimal(b)]


def interval_rows():
    mp.iv.dps = 60
    I = mp.iv
    w = [I.mpf(0)]
    for d in range(1,Q):
        r = I.mpf(min(d,Q-d))
        w.append(I.cos(I.mpf("0.18575")*r+I.mpf("0.1625"))/(1+I.power(r,I.mpf("1.8"))))
    lam = {k:sum((w[d]*(1-I.cos(2*I.pi*k*d/Q)) for d in range(1,Q)),I.mpf(0)) for k in (3,4,5,6)}
    a = [sum((lam[k]*I.cos(2*I.pi*k*d/Q)/6 for k in (3,4,5)),I.mpf(0))+lam[6]*(-1)**d/12 for d in range(Q)]
    tr = Q*a[0]
    flat = [(tr/7)*(sum((I.cos(2*I.pi*k*d/Q)/6 for k in (3,4,5)),I.mpf(0))+(-1)**d/I.mpf(12)) for d in range(Q)]
    potts = [(tr/11)*(I.mpf(int(d==0))-I.mpf(1)/12) for d in range(Q)]
    return {"FIN_A7":a,"FULL_POTTS_TRACE":potts,"FLAT_P7_TRACE":flat}


def iv_gain(lo,hi=None):
    I = mp.iv
    exact = lambda x:I.mpf(x.numerator)/I.mpf(x.denominator)
    if hi is None:
        return exact(lo)
    return I.mpf([exact(lo).a,exact(hi).b])


def iv_probability(row,g):
    I = mp.iv
    # Subtract the common diagonal log-weight: normalized law is unchanged.
    weights = [I.exp(g*(v-row[0])/2) for v in row]
    total = sum(weights,I.mpf(0))
    return [w/total for w in weights]


def tv_interval(p,q):
    lows,highs = [],[]
    for a,b in zip(p,q):
        lo,hi = endpoints(a-b)
        lows.append(max(Fraction(0),lo,-hi))
        highs.append(max(abs(lo),abs(hi)))
    return sum(lows,Fraction())/2,sum(highs,Fraction())/2


def gain_cover(rows):
    lower,step,cells = Fraction(297,100),Fraction(1,1000),60
    grid = [(lower+i*step,lower+(i+1)*step) for i in range(cells)]
    laws = {name:[iv_probability(rows[name],iv_gain(a,b)) for a,b in grid]
            for name in ("FIN_A7","FLAT_P7_TRACE")}
    public_laws = {name:[[box(z) for z in p] for p in values] for name,values in laws.items()}
    fractions = {name:[[(Fraction(lo),Fraction(hi)) for lo,hi in p] for p in values] for name,values in public_laws.items()}
    bounds = []
    for i in range(cells):
        bounds.append([sum((max(Fraction(0),a[0]-b[1],b[0]-a[1]) for a,b in zip(fractions["FIN_A7"][i],fractions["FLAT_P7_TRACE"][j])),Fraction())/2 for j in range(cells)])
    minimum = min(v for row in bounds for v in row)
    # One exact decimal witness gives an upper bound, not a minimizer proof.
    witness_g = (Fraction(297,100),Fraction(297343,100000))
    witness_laws = [[box(z) for z in iv_probability(rows[name],iv_gain(g))] for name,g in zip(("FIN_A7","FLAT_P7_TRACE"),witness_g)]
    witness_hi = sum((max(abs(Fraction(a[0])-Fraction(b[1])),abs(Fraction(a[1])-Fraction(b[0]))) for a,b in zip(*witness_laws)),Fraction())/2
    assert minimum > Fraction(1,100)
    public = Fraction(decimal(minimum,False))
    assert public > Fraction(1,100)
    return {
        "target":"inf over independent g_FIN,g_FLAT in [2.97,3.03] of ideal N2 stationary TV",
        "assumptions":["fixed exact-decimal FIN and trace-matched flat P7", "zero one-body field", "ideal difference law before detector", "correct mpmath.iv enclosures"],
        "gain_grid":{"origin":"2.97","step":"0.001","cells_per_axis":cells,"last_endpoint":"3.03","closed_boxes":cells*cells},
        "coverage":"all Cartesian products of consecutive closed gain intervals; endpoints shared, no gaps",
        "global_ideal_TV_lower":decimal(minimum,False),
        "point_witness_g":["2.97","2.97343"],
        "global_ideal_TV_upper_from_witness":decimal(witness_hi),
        "point_witness_probability_boxes":witness_laws,
        "not_concluded":"exact minimizing gain pair or exact minimax value",
        "per_hypothesis_detector_output_TV_allowance":"0.005",
        "robust_observed_law_lower_after_two_detector_allowances":decimal(public-Fraction(1,100),False),
        "robust_scope":"hypothetical model sets within output-TV .005 of each ideal law; gains may differ independently",
        "calibration_achieved":False,
        "cell_lower_bounds":[[decimal(v,False) for v in row] for row in bounds],
        "probability_boxes":[{"model":name,"boxes":values} for name,values in public_laws.items()],
    }


def quotient(row,kind,g=3):
    p = softmax(g*np.asarray(row)/2)
    B = np.zeros((Q,Q))
    for d,v in itertools.product(range(Q),repeat=2):
        if d==v:
            continue
        delta = g*(row[v]-row[d])/2
        if kind=="heat_bath": rate=2*p[v]
        elif kind=="metropolis": rate=(2/Q)*min(1,math.exp(delta))
        elif kind=="barker": rate=(4/Q)/(1+math.exp(-delta))
        else: raise ValueError(kind)
        B[d,v]=rate
    np.fill_diagonal(B,-B.sum(axis=1))
    return B


def parent(A,kind,g=3):
    L = np.zeros((Q*Q,Q*Q))
    for a,b in itertools.product(range(Q),repeat=2):
        x = a*Q+b
        for site in (0,1):
            old,other = (a,b) if site==0 else (b,a)
            fields = g*A[:,other]/2
            p = softmax(fields)
            for new in range(Q):
                if new==old:
                    continue
                dest = new*Q+b if site==0 else a*Q+new
                delta = fields[new]-fields[old]
                rate = p[new] if kind=="heat_bath" else ((1/Q)*min(1,math.exp(delta)) if kind=="metropolis" else (2/Q)/(1+math.exp(-delta)))
                L[x,dest] += rate
    np.fill_diagonal(L,-L.sum(axis=1))
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--replay",type=Path)
    args = ap.parse_args()
    files = []
    for line in (SOURCE/"MANIFEST.sha256").read_text().splitlines():
        h,name = line.split(None,1)
        assert sha(SOURCE/name)==h,name
        files.append(name)
    uncovered = sorted(str(p.relative_to(SOURCE)) for p in SOURCE.rglob("*") if p.is_file() and str(p.relative_to(SOURCE)) not in files and p.name!="MANIFEST.sha256")
    assert not uncovered,uncovered
    stage_count=0
    for stage in sorted(SOURCE.glob("RG-E0[1-5]")):
        for line in (stage/"INPUTS.sha256").read_text().splitlines():
            h,name=line.split(None,1)
            assert sha(SOURCE/name)==h,(stage,name)
            stage_count+=1
    with zipfile.ZipFile(SOURCE/"source_context/FIN_RG_D01_D05_FULL_HANDOFF_20261014.zip") as z:
        members = [n for n in z.namelist() if not n.endswith("/")]
        assert all((ROOT/n).is_file() and hashlib.sha256(z.read(n)).hexdigest()==sha(ROOT/n) for n in members)
    core = ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    spec=importlib.util.spec_from_file_location("core",core)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    A=module.build_rank7()[-1]
    all_models=module.build_countermodels(A)
    models={name:all_models[name] for name in ("FIN_A7","FULL_POTTS_TRACE","FLAT_P7_TRACE")}
    saved=read(SOURCE/"RG-E02/SHAPE_PREDICTIONS.json")
    source_dyn=read(SOURCE/"RG-E04/MULTITIME_PREDICTION_DESIGN.json")["N2_two_time"]
    d=np.array([(b-a)%Q for a,b in itertools.product(range(Q),repeat=2)])
    R=np.eye(Q)[d]
    numeric=[]
    for name,kernel in models.items():
        p=softmax(1.5*kernel[0])
        assert np.max(np.abs(p-saved["observable_laws"][name]["N2_diff_g3"]))<1e-13
        ii,jj,kk=np.indices((Q,Q,Q))
        f=(5.145228719489142/3)*(kernel[ii,jj]+kernel[ii,kk]+kernel[jj,kk])
        w=np.exp(f-f.max());vis=(w/w.sum()).sum(axis=2)
        assert np.max(np.abs(vis-saved["observable_laws"][name]["M3_ordered_piY_g5"]))<1e-13
        for kind in ("heat_bath","metropolis","barker"):
            B=quotient(kernel[0],kind)
            L=parent(kernel,kind)
            defect=float(np.max(np.abs(L@R-R@B)))
            assert defect<1e-12
            initial=np.ones(Q)/Q
            pt=initial@expm(.4*B)
            joint=pt[:,None]*expm(.6*B)
            recorded=next(r["joint_13x13"] for r in source_dyn["joint_13x13_by_model_kinetic"][name] if r["kinetic"]==kind)
            diff=float(np.max(np.abs(joint-np.asarray(recorded)[:12,:12])))
            assert diff<1e-12 and np.min(joint)>=0 and abs(joint.sum()-1)<1e-12
            assert np.max(np.abs(quotient(kernel[0],kind,0)-quotient(kernel[0],"heat_bath",0)))<1e-14
            if kind=="heat_bath":
                formula=(np.exp(-.8)*initial+(1-np.exp(-.8))*p)[:,None]*(np.exp(-1.2)*np.eye(Q)+(1-np.exp(-1.2))*p[None,:])
                assert np.max(np.abs(formula-joint))<1e-13
            numeric.append({"model":name,"kinetic":kind,"parent_quotient_intertwining_max":defect,"saved_joint_max_discrepancy":diff,"joint_mass":float(joint.sum()),"generator":B.tolist()})
    bad=A.copy();bad[0,1]+=.01;bad[1,0]+=.01
    broken=parent(bad,"heat_bath")@R
    defect=max(float(np.max(np.ptp(broken[d==v],axis=0))) for v in range(Q))
    assert defect>1e-5
    rows=interval_rows()
    cert=read(SOURCE/"RG-E03/STATIC_INTERVAL_CERTIFICATE.json")
    laws={name:iv_probability(row,iv_gain(Fraction(3))) for name,row in rows.items()}
    for name,p in laws.items():
        for z,(lo,hi) in zip(p,cert["individual_law_enclosures"][name]):
            a,b=endpoints(z)
            assert Fraction(lo)<=a<=b<=Fraction(hi)
    stat={}
    for i,name in enumerate(models):
        for other in list(models)[i+1:]:
            lo,hi=tv_interval(laws[name],laws[other])
            public=cert["pairwise_tv_public_outward"][name+"_vs_"+other]
            assert Fraction(public[0])<=lo<=hi<=Fraction(public[1])
            stat[name+"_vs_"+other]=[decimal(lo,False),decimal(hi)]
    cover=gain_cover(rows)
    args.output.with_name("GAIN_BOX_CERTIFICATE.json").write_text(json.dumps(cover,indent=2)+"\n")
    normalize=lambda s:re.sub(r"\s+","",html.unescape(s).replace("\\","").replace("*",""))
    replay=None
    if args.replay:
        replay={"directory":str(args.replay),"verifier":read(args.replay/"INDEPENDENT_VERIFICATION.json"),
                "static_certificate_byte_identical":sha(args.replay/"RG-E03/STATIC_INTERVAL_CERTIFICATE.json")==sha(SOURCE/"RG-E03/STATIC_INTERVAL_CERTIFICATE.json")}
    summary={
        "date":"2026-10-14","manifest":{"entries":len(files),"mismatches":0,"uncovered":uncovered,"stage_hashes_valid":stage_count},
        "embedded_predecessor_members_byte_matched":len(members),
        "plan_byte_identical":sha(ROOT/"FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014.md")==sha(SOURCE/"source_context/FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014_USER.md"),
        "plan_format_normalized_identical":normalize((ROOT/"FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014.md").read_text())==normalize((SOURCE/"source_context/FIN_GENERATIVE_OPERATIONAL_NEXT_CHAIN_20261014_USER.md").read_text()),
        "input_hashes":{str(core.relative_to(ROOT)):sha(core),"AGENTS.md":sha(ROOT/"AGENTS.md"),"source_manifest":sha(SOURCE/"MANIFEST.sha256"),"checker":sha(Path(__file__))},
        "numeric_quotient_checks":numeric,"symmetry_breaking_negative_control_defect":defect,
        "static_interval_rechecks":stat,
        "sharp_gain_cover_summary":{k:v for k,v in cover.items() if k not in ("cell_lower_bounds","probability_boxes")},
        "source_replay":replay,"versions":{"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__},
        "limits":["Exact quotient requires stated circulant/zero-field/equivariant kinetics", "No physical calibration/validation records", "Dynamic matrix propagation is floating numerical", "Static/cover enclosures assume mpmath.iv", "The positive robust margin is conditional on hypothetical per-law detector errors, not hardware fidelity"]}
    args.output.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({"manifest":summary["manifest"],"quotient_laws":len(numeric),"embedded_predecessor":len(members),"gain_cover":summary["sharp_gain_cover_summary"],"plan_equivalent":summary["plan_format_normalized_identical"]},indent=2))


if __name__=="__main__":
    main()
