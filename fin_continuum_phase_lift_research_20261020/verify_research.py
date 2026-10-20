#!/usr/bin/env python3
"""Exact small controls for CPL; no simulation or native count construction."""
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
    "fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md",
    "fin_reversible_bracket_averaging_review_20261016/PROOF_AUDIT_AND_SCOPE.md",
    "fin_count_partition_feasibility_review_20261016/PROOF_AUDIT_AND_CAP_ENCLOSURE.md",
    "fin_count_partition_feasibility_campaign_20261016/FP-002/POSITIVITY_AND_CAP_PROOF.md",
    "fin_stationary_kinetic_rigidity_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_auxiliary_mediator_review_20261016/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
    "fin_mesoscopic_hall_20261019/HANDOFF.md",
]
PROTECTED=[
    "fin_reversible_bracket_averaging_review_20261016",
    "fin_count_partition_feasibility_review_20261016",
    "fin_count_partition_feasibility_campaign_20261016",
    "fin_stationary_kinetic_rigidity_review_20261016",
    "fin_critical_metastability_review_20261015",
    "fin_auxiliary_mediator_review_20261016",
    "fin_mesoscopic_hall_20261019",
    "fin_rp004_full_input_self_audit_20261020",
    "fin_inverse_history_field_research_20261020",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen=set()
    for row in (folder/"MANIFEST.sha256").read_text().splitlines():
        h,name=row.split(None,1)
        name=name.removeprefix("*")
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name)==h,(folder,name)
        seen.add(name)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and p.suffix!=".pyc"}
    assert actual==seen|{"MANIFEST.sha256"},(folder,actual-seen,seen-actual)
    return {"entries":len(seen),"sha256":sha(folder/"MANIFEST.sha256")}


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))


def controls():
    import sympy as s
    checks=[]
    def ck(name,ok):
        assert name not in checks and bool(ok),name
        checks.append(name)
    h,p=s.symbols("h p",positive=True)
    z=s.symbols("z",real=True)
    a=h*(1-h)**2
    scale=s.log(h)-s.log(1-h)+1/(1-h)
    ck("scale_derivative_exact",s.simplify(s.diff(scale,h)-1/a)==0)
    ck("scale_is_harmonic",s.simplify(s.diff(a*s.diff(scale,h),h))==0)
    ck("scale_diverges_at_zero",s.limit(scale,h,0,dir="+")==-s.oo)
    ck("scale_diverges_at_one",s.limit(scale,h,1,dir="-")==s.oo)
    ck("radial_conductivity_endpoints",a.subs(h,0)==a.subs(h,1)==0 and s.diff(a,h).subs(h,0)==1 and s.diff(a,h).subs(h,1)==0)
    ck("uniform_h_generator_invariance_polynomial",s.integrate(a*2+s.diff(a,h)*2*h,(h,0,1))==0)
    RA=-2*s.log(1-p*h)
    RB=-2*s.log((1-p)*(1-h))
    jacA=(1-p*h)*s.diff(RA,h)/2
    jacB=(1-p)*(1-h)*s.diff(RB,h)/2
    ck("phase_A_Gaussian_radial_Jacobian",s.simplify(jacA-p)==0)
    ck("phase_B_Gaussian_radial_Jacobian",s.simplify(jacB-(1-p))==0)
    pp,hh=F(5,8),F(1,3)
    ua,ub=pp*hh,pp+(1-pp)*hh
    ck("phase_intrinsic_and_chart_inverse",ua<pp<ub and ua/pp==(ub-pp)/(1-pp)==hh)
    ck("mixture_radial_density_uniform",pp/pp==(1-pp)/(1-pp)==1)
    pz=s.Rational(1,2)+s.exp(-z*z)/8
    ck("original_frozen_partition_derivative",s.simplify(s.diff(pz,z)+z*s.exp(-z*z)/4)==0)
    rho=s.exp(-z**4/12)
    b=-z**3/3
    ck("quartic_zero_current_density",s.simplify(s.diff(rho,z)-b*rho)==0)
    identity=s.simplify(s.diff(z*s.exp(-z*z)*rho,z)/rho)
    ck("phase_residual_integration_identity",s.simplify(identity-s.exp(-z*z)*(1-2*z*z-z**4/3))==0)
    JA,rr=s.symbols("JA rr",positive=True)
    JB=rr-JA
    ga=2*JA/(rr*p)
    gb=2*JB/(rr*(1-p))
    ck("conditional_mean_soft_bracket_exact",s.simplify(p*ga+(1-p)*gb-2)==0)
    ck("soft_Lyapunov_uniform_bound",F(8,3)>0 and s.Poly(-2*z**4/3,z).coeff_monomial(z**4)<0)
    M=s.symbols("M",positive=True)
    eps=F(1,8)
    residualA=-4*eps*M
    residualB=4*eps*M
    ck("exact_unpostselected_residual_mass",s.simplify(residualA+residualB)==0 and residualA==-M/2 and residualB==M/2)
    ck("IBA_defect_coefficient_is_M",s.simplify(-residualA+residualB-M)==0)
    ck("coarse_exp_bound_for_positive_M",F(3)*F(12,11)==F(36,11)<4)
    ck("coarse_M_lower_bound",F(1,4)*F(1,4)/8==F(1,128))
    ck("quartic_tail_dominates_exponential_after_three",3**3>12)
    # Explicit constants follow the analytical inequalities, not a sampled domain scan.
    pmin,pmax=F(1,2),F(5,8)
    A_cap=pmax/(1-pmax)**2
    B_cap=F(2)
    ck("phase_A_radial_trace_bound",A_cap==F(40,9))
    ck("both_phase_pair_fast_trace_bound",max(A_cap,B_cap)+2==F(58,9))
    r_z_cap=F(1,32)/(2*pmin*(1-pmax)**2)
    ck("soft_radial_derivative_cap",r_z_cap==F(2,9))
    ck("soft_lift_noise_trace_cap",F(8,3)*r_z_cap==F(16,27))
    ck("whitened_total_fast_trace_cap",F(58,9)+16==F(202,9))
    ck("native_trace_growth_can_be_matched_necessarily",F(4,3)-F(5,12)*F(8,3)==F(2,9)>0)
    # All single-label increments, at a purely kinematic N144 fixture.
    N,sqrtN=144,12
    sigma=tuple(F((-1)**j) for j in range(12))
    rates={}
    G=F(0)
    fast=F(0)
    cross=same=F(0)
    edge_checks=[]
    for aa in range(12):
        for bb in range(12):
            if aa==bb:
                continue
            delta=tuple(F(j==bb)-F(j==aa) for j in range(12))
            ds=sigma[bb]-sigma[aa]
            proj=tuple(delta[j]-ds*sigma[j]/12 for j in range(12))
            rate=F((aa+2)*(bb+3),7+aa+bb)
            rates[(aa,bb)]=rate
            norm=dot(proj,proj)
            edge_checks.append(norm==(F(5,3) if ds else F(2)))
            edge_checks.append(dot(proj,sigma)==0 and sum(proj)==0)
            dz2=ds*ds/(N*sqrtN)
            dy2=norm/N
            G+=sqrtN*rate*dz2
            fast+=sqrtN*rate*dy2
            if ds:
                cross+=rate
            else:
                same+=rate
    ck("all_132_single_copy_fast_projections",len(rates)==132 and all(edge_checks))
    ck("critical_soft_trace_exact_factor_N",G==4*cross/N)
    ck("native_fast_trace_identity_exact",fast==F(5,12)*sqrtN*G+F(2,sqrtN)*same)
    Ccap=(cross+same)/N
    ck("native_trace_upper_and_lower",F(5,12)*sqrtN*G<=fast<=2*Ccap*sqrtN)
    ck("wrong_missing_sqrtN_factor_rejected",fast!=F(5,12)*G+F(2,sqrtN)*same)
    raw_cos=tuple(F((1,0,-1,0)[j%4]) for j in range(12))
    cosine_cross_checks=[]
    local_moment_slacks=[]
    for aa,bb in rates:
        ds=sigma[bb]-sigma[aa]
        local_moment_slacks.append((raw_cos[bb]-raw_cos[aa])**2-ds**2/4)
        if sigma[bb]!=sigma[aa]:
            cosine_cross_checks.append((raw_cos[bb]-raw_cos[aa])**2==1)
    ck("all_cross_parity_cosine_jumps_squared_one",len(cosine_cross_checks)==72 and all(cosine_cross_checks))
    ck("local_cosine_moment_dual_all_132_edges",min(local_moment_slacks)==0 and local_moment_slacks.count(F(0))==114 and local_moment_slacks.count(F(4))==18)
    spectral=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
    geometry=json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("native_spectral_input_matches_geometry",sha(spectral)==geometry["spectral_sha256"])
    boxes=json.loads(spectral.read_text())["exact"]["laplacian_intervals"]
    gamma_lower=1-F(boxes[3][1])/F(boxes[6][0])
    ck("actual_gamma3_exceeds_one_eighth",gamma_lower>F(1,8))
    radial_strip=F(40,9)*F(1,20)**2
    angular_strip=2*F(1,20)**2
    soft_strip=F(16,27*20)
    witness_upper=radial_strip+angular_strip+soft_strip
    ck("frozen_lift_positive_mass_cone_witness_upper",witness_upper==F(247,5400)<F(1,20))
    required_lower=F(1,8)*F(1,2)*F(8,5)
    ck("native_required_cosine_variance_lower",required_lower==F(1,10) and required_lower-witness_upper>F(1,20))
    mass_lower=F(1,2)*F(1,20)*F(1,40)
    ck("positive_witness_mass_bound",mass_lower==F(1,1600))
    ck("positive_joint_characteristic_L1_obstruction",mass_lower*F(1,20)==F(1,32000)>0)
    # A tight state process can still have exploding compensated bracket.
    n=13
    amplitude,rate=F(1,n),n**3
    ck("two_state_countercontrol_for_bracket_vs_state_tightness",amplitude**2*rate==n and amplitude<F(1,10))
    ck("non_tight_bracket_probability_bound_algebra",(F(3,2)-1)/4==F(1,8))
    f=s.exp(-z*z)
    ck("bounded_observable_second_derivative",s.simplify(s.diff(f,z,2)-(4*z*z-2)*f)==0)
    L0f=b*s.diff(f,z)+s.diff(f,z,2)
    ck("bounded_observable_generator_polynomial_times_Gaussian",s.simplify(L0f/f-(2*z**4/3+4*z*z-2))==0)
    gap=s.symbols("gap",positive=True)
    gA=2-gap
    gB=2+p*gap/(1-p)
    var=s.simplify(p*(gA-2)**2+(1-p)*(gB-2)**2)
    ck("positive_compression_defect_formula",s.simplify(var-p*gap**2/(1-p))==0)
    vfixture=var.subs({p:s.Rational(5,8),gap:s.Rational(1,3)})
    ck("exact_finite_covariance_defect_fixture",vfixture==s.Rational(5,27) and vfixture/4==s.Rational(5,108))
    # A finite variable-angle-clock chain illustrates preservation with no feedback.
    base=[[-F(1),F(1)],[F(1),-F(1)]]
    full=[[F(0)]*8 for _ in range(8)]
    for zz in range(2):
        for th in range(4):
            i=4*zz+th
            full[i][4*(1-zz)+th]+=1
            for step in (-1,1):
                full[i][4*zz+(th+step)%4]+=zz+1
            full[i][i]=-sum(full[i])
    ck("finite_variable_clock_model_reversible",all(sum(row)==0 for row in full) and all(full[i][j]==full[j][i] for i in range(8) for j in range(8)))
    ck("finite_angle_clock_has_no_base_feedback",all(sum(full[4*zz+th][4*vv+tt] for tt in range(4))==base[zz][vv] for zz in range(2) for vv in range(2) for th in range(4)))
    return checks,{"native_increment_fixture_N":N,"single_copy_edges_checked":len(rates),
                   "kinematic_fixture_scope":"increments only; no MT-domain or stationary-generator assertion",
                   "nonuniform_arbitrary_fixture_rates":True,"native_generator_constructed":False,
                   "finite_variable_clock_states":8,"finite_compression_defect":"5/108",
                   "comparison_M_strict_lower":"1/128",
                   "phase_pair_fast_trace_upper":"58/9","whitened_full_fast_trace_upper":"202/9",
                   "soft_lift_trace_upper":"16/27","trace_margin_at_kappa_sqrtN":"(2/9)sqrt(N)",
                   "native_componentwise_cone_witness":{"gamma3_lower_exceeds":"1/8","cosine_dual_edge_slacks":{"zero":114,"four":18},
                   "CPL_normalized_cosine_variance_upper_on_U":"247/5400","native_required_lower_on_U":"1/10",
                   "stationary_witness_mass_strict_lower":"1/1600","joint_characteristic_error_liminf_lower":"1/32000"}}


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
    before={n:manifest(ROOT/n) for n in PROTECTED}
    checks,details=controls()
    assert inputs=={n:sha(ROOT/n) for n in INPUTS}
    assert before=={n:manifest(ROOT/n) for n in PROTECTED}
    result={"date":"2026-10-20","status":"CONTINUUM_REALIZED_FROZEN_NATIVE_CHARACTERISTIC_TRANSFER_REJECTED",
            "model_id":"FIN-CONT-PHASE-LIFT-001; diagnostic only, not adopted FIN law",
            "bounded_checks":len(checks),"checks":checks,"finite_controls":details,
            "input_hashes":inputs,"protected_manifests":before,"originals_unchanged":True,
            "same_author":True,"independent_peer_review":False,"subagents_used":False,
            "finite_N_IBA_counterfamily":False,"native_unconditional_IBA":"OPEN",
            "frozen_CPL_native_characteristic_transfer":"REJECTED_POSITIVE_MASS_COSINE_VARIANCE_WITNESS",
            "joint_characteristic_error_liminf_lower":"1/32000",
            "general_native_phase_flow_feasibility":"OPEN",
            "F01":"OPEN_SOURCE","AGENTS_modified":False,"broad_goal":"ACTIVE_NOT_COMPLETE"}
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    links()
    input_text="".join(f"{h}  {n}\n" for n,h in inputs.items())
    if args.record:
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        files=sorted(p for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"RESULTS.json").read_text())==result
        assert (HERE/"INPUTS.sha256").read_text()==input_text
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"bounded_checks":len(checks),"inputs":len(inputs),
                      "protected_manifests":len(before),"protected_entries":sum(v["entries"] for v in before.values()),
                      "native_count_counterexample":False,"F01":"OPEN_SOURCE","same_author":True},indent=2))


if __name__=="__main__":
    main()
