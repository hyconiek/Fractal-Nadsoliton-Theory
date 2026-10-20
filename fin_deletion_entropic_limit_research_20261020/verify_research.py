#!/usr/bin/env python3
"""Bounded exact controls for DEL; not a mechanical proof of the limit."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=[
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_rbcut_weighted_degree_review_20261016/CRITICAL_DELETION_DRIFT_OBSTRUCTION.md",
    "fin_rbcut_weighted_degree_review_20261016/PROOF_AUDIT_AND_CORRECTIONS.md",
    "fin_balanced_phase_current_campaign_20261016/MARGINAL_AND_SUPPORT_CONTRACT.json",
    "fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
    "fin_stationary_kinetic_rigidity_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_stationary_kinetic_rigidity_review_20261016/REVERSIBLE_POISSON_IBA_CRITERION.md",
    "fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md",
    "fin_continuum_phase_lift_research_20261020/NATIVE_COSINE_VARIANCE_OBSTRUCTION.md",
    "fin_mesoscopic_hall_20261019/EXACT_LATTICE_COSETS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_critical_metastability_review_20261015/verification.json",
]
PROTECTED=[
    "fin_rbcut_weighted_degree_review_20261016",
    "fin_balanced_phase_current_review_20261017",
    "fin_uniform_critical_campaign_20261014",
    "fin_critical_metastability_review_20261015",
    "fin_stationary_kinetic_rigidity_review_20261016",
    "fin_reversible_bracket_averaging_review_20261016",
    "fin_continuum_phase_lift_research_20261020",
    "fin_mesoscopic_hall_20261019",
    "fin_rp004_full_input_self_audit_20261020",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen=set()
    for row in (folder/"MANIFEST.sha256").read_text().splitlines():
        digest,name=row.split(None,1)
        name=name.removeprefix("*")
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name)==digest,(folder,name)
        seen.add(name)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    assert actual==seen|{"MANIFEST.sha256"},(folder,actual-seen,seen-actual)
    return {"entries":len(seen),"sha256":sha(folder/"MANIFEST.sha256")}


def matmul(A,B):
    return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]


def matvec(A,v):
    return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]


def controls():
    import sympy as s
    checks=[]
    def ck(name,ok):
        assert name not in checks and bool(ok),name
        checks.append(name)
    spectral=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
    geometry=json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("current_spectra_match_admitted_geometry",sha(spectral)==geometry["spectral_sha256"])
    mt=json.loads((ROOT/"fin_critical_metastability_review_20261015/verification.json").read_text())
    seed="fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    ck("current_A7_seed_matches_saved_scientific_input",sha(ROOT/seed)==mt["input_hashes"][seed])
    raw=json.loads(spectral.read_text())["exact"]["laplacian_intervals"]
    ck("inherited_fast_gap_is_strict",0<F(raw[3][0])<=F(raw[3][1])<F(raw[4][0])<=F(raw[4][1])<F(raw[5][0])<=F(raw[5][1])<F(raw[6][0]))
    q=12
    sig=[F((-1)**j) for j in range(q)]
    PF=[[F(i==j)-F(1,q)-sig[i]*sig[j]/q for j in range(q)] for i in range(q)]
    Lop=[[F(6) if i==j else -F(sig[i]!=sig[j]) for j in range(q)] for i in range(q)]
    ck("opposite_graph_constant_kernel",matvec(Lop,[F(1)]*q)==[0]*q)
    ck("opposite_graph_soft_eigenvalue_12",matvec(Lop,sig)==[12*x for x in sig])
    ck("opposite_graph_fast_eigenvalue_6",matmul(Lop,PF)==[[6*x for x in row] for row in PF])
    K=[[x/F(q*q) for x in row] for row in Lop]
    ck("actual_fast_mobility_one_24",matmul(K,PF)==[[x/F(24) for x in row] for row in PF])
    ck("actual_fast_covariance_one_12",matmul([[2*x for x in row] for row in K],PF)==[[x/F(12) for x in row] for row in PF])
    ck("actual_soft_bracket_two",sum(sig[i]*x for i,x in enumerate(matvec([[2*x for x in row] for row in K],sig)))==2)
    up=[]
    down=[]
    for a in range(q):
        for b in range(q):
            if sig[a]==sig[b]:
                continue
            ds=sig[b]-sig[a]
            vec=[F(j==b)-F(j==a)-ds*sig[j]/q for j in range(q)]
            (up if ds>0 else down).append(vec)
    ck("mixed_frame_cancellation_separately_each_parity",len(up)==len(down)==36 and all(sum(v[j] for v in up)==sum(v[j] for v in down)==0 for j in range(q)))
    m,t=s.symbols("m t",real=True)
    for sign in (-1,1):
        bj=(1-sign*m)/2*(1+sign*t)/q-(1+sign*m)/q*(1-sign*t)/2
        ck(f"exact_meanfield_parity_line_sign_{sign}",s.expand(bj-sign*(t-m)/q)==0)
    w,Aw,g=s.symbols("w Aw g",real=True)
    ck("fast_linearization_factor_24",s.expand((g*Aw/q-w)/2+(q*w-g*Aw)/(2*q))==0)
    ck("critical_cubic_not_refitted",s.series(s.tanh(m)-m,m,0,5).removeO()==-m**3/3)
    ck("soft_third_characteristic_remainder_exponent",F(3,2)-3*F(3,4)==-F(3,4))
    ck("soft_dual_calibration_remainder_exponent",F(3,2)-F(3,4)-2*F(1,2)==-F(1,4))
    ck("fast_dual_remainder_exponent",F(3,2)-3*F(1,2)-F(1,2)==-F(1,2))
    ck("fast_Cauchy_collapse_exponent",F(1,4)-F(1,2)==-F(1,4))
    yy=s.symbols("yy",real=True)
    harmonic_derivative=s.exp(yy*yy/2)
    ck("fast_harmonic_alone_need_not_be_constant",s.simplify(s.diff(harmonic_derivative,yy)-yy*harmonic_derivative)==0 and harmonic_derivative.subs(yy,0)==1)
    ck("full_scaled_cell_volume_product",6*2==12 and F(5)+F(3,4)==F(23,4))
    z=s.symbols("z",real=True)
    rho=s.exp(-z**4/12)
    b=-z**3/3
    p=s.Rational(1,2)+s.exp(-z*z)/8
    pB=1-p
    dA=s.diff(p,z)/p
    dB=s.diff(pB,z)/pB
    bA,bB=b+dA,b+dB
    ck("phase_A_zero_current_limit",s.simplify(s.diff(rho*p,z)-bA*rho*p)==0)
    ck("phase_B_zero_current_limit",s.simplify(s.diff(rho*pB,z)-bB*rho*pB)==0)
    ck("global_conditional_drift_cancels",s.simplify(p*bA+pB*bB-b)==0)
    ck("phase_A_linear_force_at_zero",s.diff(bA,z).subs(z,0)==-s.Rational(2,5))
    ck("phase_B_linear_force_at_zero",s.diff(bB,z).subs(z,0)==s.Rational(2,3))
    ck("old_weighted_boundary_signs_recovered",s.simplify(z*rho*s.diff(p,z)+z*z*rho*s.exp(-z*z)/4)==0)
    # Exact edge-field duality on a separate small reversible weighted graph.
    pi=[F(1,6),F(1,3),F(1,2)]
    c=[[F(0),F(1,7),F(1,11)],[F(1,7),F(0),F(2,5)],[F(1,11),F(2,5),F(0)]]
    vals=[F(1),F(-2),F(3)]
    v=[[F(0),F(2,3),F(5,8)],[-F(2,3),F(0),-F(1,4)],[-F(5,8),F(1,4),F(0)]]
    E=sum(c[i][j]*(vals[j]-vals[i])**2 for i in range(3) for j in range(3))/2
    B=sum(c[i][j]*(vals[j]-vals[i])*v[i][j] for i in range(3) for j in range(3))/2
    Q=sum(c[i][j]*v[i][j]**2 for i in range(3) for j in range(3))/2
    divergence=-sum(pi[i]*vals[i]*sum(c[i][j]/pi[i]*v[i][j] for j in range(3)) for i in range(3))
    ck("exact_antisymmetric_edge_divergence_sign",B==divergence)
    square=sum(c[i][j]*((vals[j]-vals[i])-v[i][j])**2 for i in range(3) for j in range(3))/2
    ck("global_energy_calibration_complete_square",E-2*B+Q==square>=0)
    ck("wrong_divergence_sign_rejected",B!=-divergence)
    # Exact cell isometry algebra; this is not an actual MT numerical fixture.
    weights=[F(2,7),F(5,7)]
    reference=[F(1,4),F(3,4)]
    f=[F(2),F(-3)]
    plain=sum(reference[i]*f[i]**2 for i in range(2))
    mapped=sum(reference[i]*(weights[i]/reference[i])*f[i]**2 for i in range(2))
    target=sum(weights[i]*f[i]**2 for i in range(2))
    ck("cell_density_multiplier_isometry",mapped==target and plain!=target)
    # Nondense slow range: a fast eigenmode resolvent collapses to zero.
    r=F(2)
    ck("nondense_range_resolvent_not_full_identity",F(1,r+2*100)<F(1,100) and F(1,r)>0)
    dF,A=s.symbols("dF A",real=True)
    M=dF-A
    reverse_increment=-dF-A
    ck("forward_backward_state_increment_MINUS",s.expand((M-reverse_increment)/2-dF)==0)
    ck("PLUS_is_not_state_increment",s.expand((M+reverse_increment)/2+A)==0)
    # A finite phase/volume model checks memory despite identical compression.
    # states 0,1,2 in A; 3,4,5 in B; z=(0,0,1,0,1,1).
    zz=[F(0),F(0),F(1),F(0),F(1),F(1)]
    phase=[0,0,0,1,1,1]
    L=[[F(0)]*6 for _ in range(6)]
    for i,j in ((0,2),(1,2),(3,4),(3,5),(0,1),(4,5)):
        L[i][j]=L[j][i]=1
    for i in range(6):
        L[i][i]=-sum(L[i])
    drift=matvec(L,zz)
    ck("finite_same_phase_volume_drift",drift==[1,1,-2,2,-1,-1])
    ck("finite_reversible_comparison",all(L[i][j]==L[j][i] for i in range(6) for j in range(6)) and all(sum(row)==0 for row in L))
    average_checks=[]
    for softlevel in (0,1):
        ids=[i for i in range(6) if zz[i]==softlevel]
        av=sum(drift[i] for i in ids)/len(ids)
        average_checks.append(av==(F(4,3) if softlevel==0 else -F(4,3)))
    ck("finite_conditional_mean_keeps_compressed_drift",all(average_checks))
    centered=[x-F(1,2) for x in zz]
    Lf=matvec(L,centered)
    L0f=[(F(4,3) if x==0 else -F(4,3)) for x in zz]
    defect=sum(x*x for x in Lf)/6-sum(x*x for x in L0f)/6
    ck("finite_second_order_memory_defect",defect==F(2,9)>0)
    # Add a reversible phase leak at an explicitly declared endpoint.
    leaking=[row.copy() for row in L]
    eps=F(1,101)
    leaking[2][3]+=eps
    leaking[3][2]+=eps
    leaking[2][2]-=eps
    leaking[3][3]-=eps
    ind=[F(x==0) for x in phase]
    test=[F(2),F(-1),F(3),F(5),F(7),F(-2)]
    a1=sum(ind[i]*x for i,x in enumerate(matvec(leaking,test)))/6
    a2=sum(test[i]*x for i,x in enumerate(matvec(leaking,ind)))/6
    ck("leaking_phase_selfadjoint_stationarity_identity",a1==a2)
    # Actual limit-law h transform, using a general positive p rather than fitting it.
    pf=s.Function("pf")(z)
    ff=s.Function("ff")(z)
    hf=s.sqrt(pf)
    L0=lambda u:s.diff(u,z,2)+b*s.diff(u,z)
    conjugate=s.simplify(L0(hf*ff)/hf-ff*L0(hf)/hf)
    expected=s.diff(ff,z,2)+(b+s.diff(pf,z)/pf)*s.diff(ff,z)
    ck("continuous_likelihood_h_transform_factor_two",s.simplify(conjugate-expected)==0)
    pp,dp=s.symbols("pp dp",positive=True)
    info=s.simplify(pp*(dp/pp)**2+(1-pp)*(-dp/(1-pp))**2)
    ck("limit_likelihood_phase_Fisher_identity",s.simplify(info-dp*dp/(pp*(1-pp)))==0)
    ck("limit_KL_bound_arithmetic",F(1,32)/F(15,64)/4==F(1,30))
    return checks,{"complete_bipartite_vertices":q,"fast_mobility":"1/24","fast_covariance":"1/12",
                   "soft_bracket":"2","cell_volume":"12N^-23/4",
                   "finite_memory_fixture_states":6,"finite_memory_defect":"2/9",
                   "finite_fixtures_are_not_MT_domain_or_native_sequence":True,
                   "limit_record_KL_upper":"S/30; limit diffusions only",
                   "full_count_graph_enumerated":False,"independent_review":False}


def links():
    for p in HERE.glob("*.md"):
        out=[]
        fence=False
        for line in p.read_text().splitlines():
            assert line==line.rstrip(),p.name
            if line.lstrip().startswith("```"):
                fence=not fence
            elif not fence and not line.startswith(("    ","\t")):
                out.append(line)
        txt=re.sub(r"`[^`\n]*`","","\n".join(out))
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)",txt):
            target=target.strip("<>").split("#")[0]
            if target and not target.startswith(("http:","https:")):
                assert (p.parent/target).exists(),(p.name,target)


def main():
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record",action="store_true")
    mode.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    if args.check:
        manifest(HERE)
    inputs={n:sha(ROOT/n) for n in INPUTS}
    protected={n:manifest(ROOT/n) for n in PROTECTED}
    checks,details=controls()
    assert inputs=={n:sha(ROOT/n) for n in INPUTS}
    assert protected=={n:manifest(ROOT/n) for n in PROTECTED}
    result={"date":"2026-10-20","status":"CONDITIONAL_STATIONARY_DELETION_SOFT_PHASE_PROCESS_LIMIT",
            "bounded_checks":len(checks),"checks":checks,"finite_controls":details,
            "input_hashes":inputs,"protected_manifests":protected,"originals_unchanged":True,
            "same_author":True,"independent_peer_review":False,"subagents_used":False,
            "primary_HB_replaced":False,"exact_KR_drift_matching":False,
            "native_unconditional_IBA":"OPEN","F01":"OPEN_SOURCE","AGENTS_modified":False,
            "broad_goal":"ACTIVE_NOT_COMPLETE"}
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    links()
    intext="".join(f"{h}  {n}\n" for n,h in inputs.items())
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        files=sorted(p for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"RESULTS.json").read_text())==result
        assert (HERE/"INPUTS.sha256").read_text()==intext
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"bounded_checks":len(checks),"inputs":len(inputs),
                      "protected_manifests":len(protected),"protected_entries":sum(v["entries"] for v in protected.values()),
                      "KR_counterfamily":False,"F01":"OPEN_SOURCE","same_author":True},indent=2))


if __name__=="__main__":
    main()
