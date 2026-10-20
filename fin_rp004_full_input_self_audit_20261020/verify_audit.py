#!/usr/bin/env python3
"""Bounded full-input self-audit controls; no independent PDE certification."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/"fin_rp004_flux_continuation_20261019"
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "FIN_RESEARCH_STRATEGY.md", "FIN_RESEARCH_EXECUTION_PLAN.md",
    "FIN_POST_RP004_RESEARCH_DECISION_20261019.md",
    "FIN_POST_RP004_EXECUTION_STATUS_20261019.json",
    "research_tasks/RP004_SHARED_CONTRACT.json",
    "fin_rp004_flux_review_20261019/MOVING_CONE_THEOREM.md",
    "fin_rp004_flux_review_20261019/BOOTSTRAP_COMPLETION.md",
    "fin_rp004_flux_review_20261019/OUTWARD_FLUX_THEOREM.md",
    "fin_rp004_flux_review_20261019/COMPLETE_FLUX_BUDGET.json",
    "fin_rp004_flux_review_20261019/INTAKE_20261019.md",
    "fin_rp004_flux_review_20261019/verify_intake.py",
    "fin_rp004_flux_review_20261019/finalize_review.py",
    "fin_rp004_independent_audit_20261019/PROOF_OR_OPEN_LEMMA.md",
    "fin_rp004_independent_audit_20261019/EXACT_GLUING_LEMMA.md",
    "fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_coherent_mode_radiation_review_20261019/RETARDED_FLUX_PROOF_ROUTE.md",
    "fin_nonlinear_stability_transfer_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_agents_causal_update_20261020/verify_update.py",
    "fin_agents_causal_update_20261020/UPDATE_CONTRACT.json",
    "fin_global_activity_defect_self_stress_20261020/USER_DIRECTION.md",
]
PROTECTED = [
    "fin_rp004_flux_review_20261019", "fin_rp004_independent_audit_20261019",
    "fin_coherent_mode_radiation_review_20261019",
    "fin_nonlinear_stability_transfer_review_20261017",
    "fin_intrinsic_entropy_source_review_20261019",
    "fin_agents_causal_update_20261020",
    "fin_finite_history_source_research_20261020",
    "fin_inverse_history_field_research_20261020",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(folder):
    return {p.relative_to(folder).as_posix():sha(p) for p in sorted(folder.rglob("*"))
            if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}


def manifest(folder):
    seen=set()
    for row in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=row.split("  ",1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name)==h,(folder,name)
        seen.add(name)
    actual=set(snapshot(folder))
    assert actual==seen|{"MANIFEST.sha256"},folder
    return {"entries":len(seen),"sha256":sha(folder/"MANIFEST.sha256")}


def controls():
    import sympy as s
    checks=[]
    def ck(name,ok):
        assert name not in checks and bool(ok),name
        checks.append(name)
    y,k=s.symbols("y k",real=True)
    h=s.Function("h")(y)
    dx=lambda f:(1-y*y)*s.diff(f,y)
    free=lambda f:-dx(dx(f))
    D=lambda n,f:dx(f)+n*y*f
    Ds=lambda n,f:-dx(f)+n*y*f
    L=lambda f:free(f)+(4-6*(1-y*y))*f
    T=lambda f:Ds(2,Ds(1,f))
    Tad=lambda f:D(1,D(2,f))
    zero=lambda f:s.simplify(s.expand(f))==0
    ck("Darboux_L_factorization",zero(Ds(2,D(2,h))-L(h)))
    ck("Darboux_intertwining",zero(L(T(h))-T(free(h)+4*h)))
    ck("Darboux_normalization_operator",zero(Tad(T(h))-(free(free(h))+5*free(h)+4*h)))
    ck("translation_in_missing_kernel",zero(Tad(1-y*y)))
    ck("shape_in_missing_kernel",zero(Tad(y*s.sqrt(1-y*y))))
    P=3*y*y-3*s.I*k*y-(k*k+1)
    dx_plane=lambda f:dx(f)+s.I*k*f
    ck("Jost_polynomial_ODE",zero(-dx_plane(dx_plane(P))+(4-6*(1-y*y))*P-(k*k+4)*P))
    ck("positive_infinity_Jost_phase",zero(P.subs(y,1)-(1-s.I*k)*(2-s.I*k)))
    F4=s.symbols("F4",positive=True)
    F6=(k*k+16)*F4/20
    transformed=3*(-s.I*k*F4/4+s.I*k*F6/6)+3*s.I*k*(F4-F6)-(k*k+1)*(-s.I*k*F4/4)
    ck("conjugated_source_transform_sign",zero(transformed-s.I*k*(k*k+4)*F4/8))
    p,nu,Om=s.symbols("p nu Om",positive=True)
    laplace=-(1/(p+s.I*Om)-p/(p*p+nu*nu)+s.I*Om/(p*p+nu*nu))/(nu*nu-Om*Om)
    ck("retarded_oscillator_Laplace_transform",zero(laplace+1/((p+s.I*Om)*(p*p+nu*nu))))
    q=s.symbols("q",positive=True)
    ck("both_root_partial_fraction",zero(1/(k*k-q*q)-(1/(k-q)-1/(k+q))/(2*q)))
    ck("outgoing_Green_jump",s.simplify(2*s.I*q*s.I/(2*q))==-1)
    ck("opposite_Green_sign_rejected",s.simplify(2*s.I*q*(-s.I)/(2*q))!=-1)
    ck("explicit_pole_cutoff_margin",(15*12)**2>25*1249)
    ck("stationary_support_disjoint_from_poles",F(1,2)<F(3,4)**2 and F(3,4)<F(15,16))
    ck("phase_curvature_positive_compact_floor",F(144)>F(125))
    # Exact smooth-profile polynomial integrals, not numerical quadrature.
    f2=s.Rational(3,4)*(1-y*y)
    norm0=s.integrate(f2,(y,-1,1))
    norm1=s.integrate(4*y*y*f2,(y,-1,1))
    norm2=s.integrate((-2+6*y*y)**2*f2,(y,-1,1))
    ck("prepared_bound_mode_normalized",norm0==1)
    ck("prepared_mode_first_derivative_norm",norm1==s.Rational(4,5))
    ck("prepared_mode_second_derivative_norm",norm2==s.Rational(16,7))
    Cv2=norm0+norm1+norm2+s.Rational(7,4)*(norm0+norm1)
    ck("uniform_Cv_three",Cv2==s.Rational(1013,140)<9)
    S_H1=s.integrate(s.Rational(81,16)*(1-y*y)**3*(y*y+(1-5*y*y)**2),(y,-1,1))
    ck("exact_source_H1_bound",0<S_H1<100)
    ck("strong_position_constant_seven",F(16,25)+F(116,25)+F(841,25)==F(973,25)<49)
    ck("strong_velocity_constant_three",F(4,5)+F(29,5)==F(33,5)<9)
    ck("Duhamel_constant_twenty_one",7*3==21 and F(6)<9)
    C1a=12*2*3*421*4
    C1b=6*3**2*421**2*4
    ck("Fr_first_coefficient",C1a==121248)
    ck("Fr_second_coefficient",C1b==38284056)
    ck("Fr_signed_self_coefficients",3**2*4**3==576 and 3**4*4**5==82944)
    pure_w=6*3**2*421*4**2+6*2*3*421**2+2*3**2*421**3
    ck("Fw_pure_coefficient_exact",pure_w==1349876718)
    ck("Fw_uniform_C2",max(12*2*3*3,6*2*3,pure_w)<1400000000)
    r,w=s.symbols("r w",positive=True)
    den=1+4*r+4*w+r*r
    for label,target in (("A",4*r),("B",4*w),("Cr_squared",r*r)):
        ck("common_K0_pays_"+label,all(v>0 for v in s.Poly(s.expand(den-target),r,w).coeffs()))
    ep,Kx,ut,ux,wt,wx,vt,vx,rt,rx=s.symbols("ep Kx ut ux wt wx vt vx rt rx")
    current=-(ep*ep*ut+wt)*(Kx+ep*ep*ux+wx)-(ep*vt+rt)*(ep*vx+rx)
    inventory=-ep**2*Kx*ut-ep**4*ut*ux-wt*Kx-ep**2*wt*ux-ep**2*ut*wx-wt*wx-ep**2*vt*vx-ep*vt*rx-ep*rt*vx-rt*rx
    ck("complete_original_current_ten_terms",zero(current-inventory) and len(s.Poly(s.expand(current),ep,Kx,ut,ux,wt,wx,vt,vx,rt,rx).terms())==10)
    ck("bound_current_omission_rejected",not zero(current-(inventory+ep**2*vt*vx)))
    ck("D1_both_endpoints",2*4*420/F(2)==1680)
    ck("D2_both_endpoints",2*6*8==96)
    ck("D3_tail_terms_dominated",2*(6+8)/F(3)<10)
    ck("D4_both_cross_terms",2*2*420==1680)
    # Differential energy identity independently checks original work signs.
    ph,ps,a,g4,g6=s.symbols("ph ps a g4 g6",real=True)
    pt,px,st,sx,pxx,sxx,pxt,sxt=s.symbols("pt px st sx pxx sxx pxt sxt",real=True)
    V=(ph*ph-1)**2/2+(a+6*ph*ph)*ps*ps/2+g4*ps**4/4+g6*ps**6/6
    et=pt*(pxx-s.diff(V,ph))+px*pxt+st*(sxx-s.diff(V,ps))+sx*sxt+s.diff(V,ph)*pt+s.diff(V,ps)*st
    jx=-pxt*px-pt*pxx-sxt*sx-st*sxx
    ck("original_local_energy_has_no_external_work",zero(et+jx))
    ck("moving_boundary_term_cannot_be_omitted",-F(1,3)+F(1,6)*2==0 and -F(1,3)!=0)
    Ls=s.symbols("Ls",positive=True)
    tau=1+12*Ls
    ck("bootstrap_log_derivative_upper",zero(3/Ls-36/tau-3/(Ls*tau)))
    ck("first_error_log_derivative_upper",zero(4/Ls-(60/tau-1/Ls)-5/(Ls*tau)))
    ck("second_error_log_derivative_upper",zero(6/Ls-(84/tau-1/Ls)-7/(Ls*tau)))
    ck("third_error_log_derivative_negative",zero(12/tau-1/Ls+1/(Ls*tau)))
    ck("whole_threshold_monotonic_after_three",-2+F(3,3)<0 and -2+F(4,3)<0 and -4+F(6,3)<0)
    c0=F(5,21)*F(125,73728)
    ck("reference_coefficient_exact",c0==F(625,1548288))
    ck("original_final_flux_coefficient_exact",c0/8==F(625,12386304))
    # A directly testable bootstrap fixture does not determine physical g4,g6.
    cr=1+21*(38405304+576*2+82944*3)
    cw=1+21*1400000000*(cr+2)
    K0=F(1,16*(1+4*cr+4*cw+cr*cr))
    ck("finite_bootstrap_fixture_no_circularity",4*cr*K0<F(1,16) and 4*cw*K0<F(1,16) and cr*cr*K0<F(1,16))
    return checks,{"S_H1_squared":str(S_H1),"v_graph_norm_squared_upper":str(Cv2),
                   "linear_Duhamel_constant":21,"Cu":420,"Cv":3,
                   "Fr_coefficient":"38405304+576M4+82944M6","Fw_coefficient":1400000000,
                   "C0":str(c0),"final_flux_coefficient":str(c0/8),
                   "fixture_scope":"M4=2,M6=3 arithmetic only, not fitted physical coefficients",
                   "no_numeric_uniform_cone_or_epsilon_threshold_claim":True}


def linkcheck():
    for p in HERE.glob("*.md"):
        prose=[]
        fence=False
        for line in p.read_text().splitlines():
            assert line==line.rstrip(),p.name
            if line.lstrip().startswith("```"):
                fence=not fence
            elif not fence and not line.startswith(("    ","\t")):
                prose.append(line)
        text=re.sub(r"`[^`\n]*`","","\n".join(prose))
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)",text):
            target=target.strip("<>").split("#")[0]
            if target and not target.startswith(("http:","https:")):
                assert (p.parent/target).exists(),(p.name,target)


def source_logs(source,context):
    assert source["status"]=="PASS_BOUNDED_CONTROLS_WITH_DISCLOSED_SOURCE_OMISSION"
    assert source["source_preserved"]==snapshot(SOURCE)
    assert source["live_AGENTS_sha256"]==sha(ROOT/"AGENTS.md")
    assert source["manifest_checks"]["."]["matching"]==10
    assert source["manifest_checks"]["."]["missing"]==["finalize.py"]
    assert source["intake_bounded_check_count"]==len(source["intake_bounded_checks"])==26
    assert len(source["copy_only_source_replay"])==2
    assert all(row["returncode"]==0 for row in source["copy_only_source_replay"])
    old=source["original_inputs"]["."]
    assert len(old)==10 and sum(row["state"]=="MATCH" for row in old)==7
    assert sum(row["state"]=="DIFFERENT_LOCAL_COUNTERPART" for row in old)==3
    update=json.loads((ROOT/"fin_agents_causal_update_20261020/UPDATE_CONTRACT.json").read_text())
    data=(ROOT/"AGENTS.md").read_bytes()
    assert hashlib.sha256(data[:update["historical_bytes"]]).hexdigest()==update["historical_sha256"]
    old_agents=next(row for row in old if row["local_path"]=="AGENTS.md")
    assert old_agents["original_sha256"]==update["historical_sha256"]
    assert context["status"]=="PASS_SCOPED_APPEND_ONLY_INSTRUCTION_UPDATE"
    logs=context["archived_context_only_replays"]
    assert len(logs)==7 and all(v["returncode"]==0 for v in logs.values())
    assert logs["fin_rp004_flux_review_20261019"]["returncode"]==0
    return {"present_original_source_files":len(source["source_preserved"]),
            "root_manifest_matches":10,"root_manifest_missing":["finalize.py"],
            "task_manifest_matches":20,"current_original_full_matches":7,
            "different_current_original_counterparts":3,"producer_assert_statements":source["producer_symbolic_assert_statements"],
            "original_AGENTS_matches_verified_historical_prefix_not_current_full_file":True,
            "source_scripts_replayed":2,"existing_review_bounded_controls":26,
            "existing_finite_cone_diagnostics":len(source["new_finite_time_spectral_diagnostic"]["rows"]),
            "historical_context_package_check":"PASS_ARCHIVED_CONTEXT_ONLY",
            "independent_review":False}


def main():
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record",action="store_true")
    mode.add_argument("--check",action="store_true")
    ap.add_argument("--source-replay",type=Path)
    ap.add_argument("--context-replay",type=Path)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    if args.check:
        manifest(HERE)
    inputs={n:sha(ROOT/n) for n in INPUTS}
    before={n:manifest(ROOT/n) for n in PROTECTED}
    source=json.loads((args.source_replay or HERE/"SOURCE_REPLAY.json").read_text())
    context=json.loads((args.context_replay or HERE/"ARCHIVED_CONTEXT_REPLAY.json").read_text())
    replay=source_logs(source,context)
    checks,constants=controls()
    result={"date":"2026-10-20","status":"ADMIT_CONDITIONAL_FINITE_WINDOW_THEOREM_SAME_AUTHOR_FULL_INPUT_AUDIT",
            "bounded_checks":len(checks),"checks":checks,"constant_payments":constants,
            "replay_scope":replay,"input_hashes":inputs,"protected_manifests":before,
            "source_preservation":True,"same_author":True,"independent_peer_review":False,
            "physical_validation":False,"F01":"OPEN_SOURCE","AGENTS_modified":False,
            "old_missing_materials_audit_modified":False,"broad_goal":"ACTIVE_NOT_COMPLETE"}
    assert inputs=={n:sha(ROOT/n) for n in INPUTS}
    assert before=={n:manifest(ROOT/n) for n in PROTECTED}
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    intext="".join(f"{h}  {n}\n" for n,h in inputs.items())
    if args.record:
        assert args.source_replay and args.context_replay
        (HERE/"SOURCE_REPLAY.json").write_text(json.dumps(source,indent=2)+"\n")
        (HERE/"ARCHIVED_CONTEXT_REPLAY.json").write_text(json.dumps(context,indent=2)+"\n")
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(intext)
        linkcheck()
        files=sorted(p for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert (HERE/"INPUTS.sha256").read_text()==intext
        assert json.loads((HERE/"RESULTS.json").read_text())==result
        linkcheck()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"bounded_checks":len(checks),"inputs":len(inputs),
                      "protected_manifests":len(before),"protected_entries":sum(v["entries"] for v in before.values()),
                      "source_root":"10 match plus absent finalize.py","source_science_replay":"2 scripts; temporary copy",
                      "historical_package_check":"ARCHIVED_CONTEXT_ONLY","independent_review":False},indent=2))


if __name__=="__main__":
    main()
