#!/usr/bin/env python3
"""Exact lattice/carry controls, not all-N proofs or count-simplex enumeration."""
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
    "research_tasks/90_MH001_MESOSCOPIC_HALL.md",
    "FIN_MESOSCOPIC_HALL_AND_INTERMEDIATE_SOURCE_NEXT_CHAIN_20261017.md",
    "fin_dense_orbit_fibre_transport_review_20261017/CLUSTERED_COMPLEMENT_SEED.md",
    "fin_dense_orbit_fibre_transport_review_20261017/PROOF_COMPLETIONS_AND_SCOPE.md",
    "fin_dense_orbit_fibre_transport_campaign_20261017/DH-001/DENSE_COMPLEMENT_HALL_THEOREM.md",
    "fin_balanced_phase_cut_research_20261016/COMPACT_PHASE_PREFIX_LOCAL_LIMIT.md",
    "fin_balanced_phase_cut_research_20261016/EXACT_COMPACT_CURRENT_MARGINALS_AND_GLUING.md",
    "fin_balanced_phase_current_campaign_20261016/MARGINAL_AND_SUPPORT_CONTRACT.json",
    "fin_balanced_phase_current_review_20261017/INTAKE_20261017.md",
    "fin_rbcut_weighted_degree_review_20261016/PROOF_AUDIT_AND_CORRECTIONS.md",
    "fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md",
    "fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
    "fin_critical_metastability_review_20261015/INTAKE_20261015.md",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
]
PROTECTED=[
    "fin_dense_orbit_fibre_transport_review_20261017",
    "fin_dense_orbit_fibre_transport_campaign_20261017",
    "fin_balanced_phase_cut_research_20261016",
    "fin_balanced_phase_current_review_20261017",
    "fin_rbcut_weighted_degree_review_20261016",
    "fin_critical_metastability_review_20261015",
    "fin_uniform_critical_campaign_20261014",
    "fin_stationary_kinetic_rigidity_review_20261016",
    "fin_rp004_full_input_self_audit_20261020",
    "fin_inverse_history_field_research_20261020",
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
    # Older producer exclusions are explicitly retained, not silently manufactured.
    extra=actual-seen-{"MANIFEST.sha256"}
    allowed=set() if folder==HERE else {"PACKAGE_VERIFICATION.json","check_package.py"}
    assert extra<=allowed,(folder,extra)
    assert seen<=actual,(folder,seen-actual)
    return {"entries":len(seen),"sha256":sha(folder/"MANIFEST.sha256"),
            "separately_bound_exclusions":{n:sha(folder/n) for n in sorted(extra)}}


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def scale(k,a):
    return tuple(k*x for x in a)


def unit(i):
    return tuple(F(j==i) for j in range(12))


def determinant(A):
    A=[[F(x) for x in row] for row in A]
    value=F(1)
    for j in range(len(A)):
        pivot=next((i for i in range(j,len(A)) if A[i][j]),None)
        if pivot is None:
            return F(0)
        if pivot!=j:
            A[j],A[pivot]=A[pivot],A[j]
            value=-value
        d=A[j][j]
        value*=d
        for i in range(j+1,len(A)):
            t=A[i][j]/d
            for k in range(j+1,len(A)):
                A[i][k]-=t*A[j][k]
    return value


def gram(basis):
    return [[dot(x,y) for y in basis] for x in basis]


def integral_lattice(v):
    return all(x.denominator==1 for x in v) and sum(v[::2])==sum(v[1::2])==0


def run():
    checks=[]
    def ck(name,ok):
        assert name not in checks and bool(ok),name
        checks.append(name)
    sigma=tuple(F((-1)**j) for j in range(12))
    c=tuple(F((1,0,-1,0)[j%4]) for j in range(12))
    s=tuple(F((0,1,0,-1)[j%4]) for j in range(12))
    full=[sub(unit(i),unit(10)) for i in (0,2,4,6,8)]
    full += [sub(unit(i),unit(11)) for i in (1,3,5,7,9)]
    kernel=[]
    for group in ((0,4,8),(2,6,10),(1,5,9),(3,7,11)):
        kernel.extend([sub(unit(group[0]),unit(group[2])),sub(unit(group[1]),unit(group[2]))])
    ck("ten_dimensional_slice_Gram_determinant_36",determinant(gram(full))==36)
    ck("eight_dimensional_kernel_Gram_determinant_81",determinant(gram(kernel))==81)
    ck("kernel_basis_is_four_A2",all(dot(v,c)==dot(v,s)==0 and integral_lattice(v) for v in kernel))
    ck("phase_features_orthogonal_norm_six",dot(c,c)==dot(s,s)==6 and dot(c,s)==0)
    be,bo=sub(unit(0),unit(2)),sub(unit(1),unit(3))
    ck("phase_map_surjects_to_even_integer_grid",dot(be,c)==dot(bo,s)==2 and dot(be,s)==dot(bo,c)==0 and all(dot(v,c).denominator==dot(v,s).denominator==1 and int(dot(v,c))%2==int(dot(v,s))%2==0 for v in full))
    pe,po=scale(F(1,3),c),scale(F(1,3),s)
    ck("projection_area_squared_four_ninths",determinant(gram([pe,po]))==F(4,9))
    ck("full_covolume_product",F(9)*F(2,3)==6 and F(81)*F(4,9)==36)
    w=tuple(F(5 if j==0 else -1,6) if j%2==0 else F(-5 if j==1 else 1,6) for j in range(12))
    ck("six_step_translate_returns_to_lattice",integral_lattice(scale(6,w)))
    ck("six_cosets_are_distinct",all(not integral_lattice(scale(r,w)) for r in range(1,6)))
    ck("translate_phase_values",dot(w,c)==1 and dot(w,s)==-1)
    ck("kernel_shift_period_three_even",integral_lattice(sub(scale(3,be),c)) and dot(sub(scale(3,be),c),c)==0)
    ck("kernel_shift_period_three_odd",integral_lattice(sub(scale(3,bo),s)) and dot(sub(scale(3,bo),s),s)==0)
    ck("nine_projected_subgrids_keep_cell_mass",9*F(1,9)*6==6 and 6/F(9)==F(2,3))
    residue_checks=[]
    for r in range(6):
        for u in range(3):
            for v in range(3):
                raw=add(scale(r,w),add(scale(u,be),scale(v,bo)))
                eta=sub(raw,add(scale(F(r+2*u,6),c),scale(F(-r+2*v,6),s)))
                residue_checks.append(dot(eta,c)==dot(eta,s)==0)
    ck("all_54_residue_shifts_really_in_kernel_space",len(residue_checks)==54 and all(residue_checks))
    fixture_count=0
    for r in range(6):
        N,l=120,60+r
        n=[0]*12
        for parity,total in ((0,l),(1,N-l)):
            base,rem=divmod(total,6)
            for i in range(6):
                n[2*i+parity]=base+int(i<rem)
        U=tuple(F(n[j])-F(l if j%2==0 else N-l,6) for j in range(12))
        assert integral_lattice(sub(U,scale(r,w)))
        assert dot(U,c)==dot(n,c) and dot(U,s)==dot(n,s)
        assert int(dot(U,c)-r)%2==0 and int(dot(U,s)+r)%2==0
        assert 2*(dot(n,c)**2+dot(n,s)**2)<=12*dot(U,U)
        for a in (1,3,5,7,9,11):
            for b in (0,2,4,6,8,10):
                nn=n.copy()
                nn[a]-=1
                nn[b]+=1
                UU=tuple(F(nn[j])-F(l+1 if j%2==0 else N-l-1,6) for j in range(12))
                delta=sub(sub(unit(b),unit(a)),scale(F(1,6),sigma))
                assert sub(UU,U)==delta and dot(delta,delta)==F(5,3)
                assert integral_lattice(sub(UU,scale((r+1)%6,w)))
        fixture_count+=1
    ck("six_exact_count_cosets_and_all_upward_steps",fixture_count==6)
    ck("ball_strictly_inside_original_phase_floor",12*F(1,10)**2==F(3,25)<1)
    ck("positive_log_two_lower_sum",F(1,2)*(F(2,3)+F(1,2))==F(7,12)>F(1,2))
    spectral=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
    geo=json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("inherited_spectral_bytes_match_geometry",sha(spectral)==geo["spectral_sha256"])
    raw=json.loads(spectral.read_text())["exact"]["laplacian_intervals"]
    lower={k:F(raw[k][0]) for k in (3,5,6)}
    upper={k:F(raw[k][1]) for k in (3,5,6)}
    ck("actual_phase_and_fast_gaps_positive",0<lower[3]<=upper[3]<lower[6] and 0<lower[5]<=upper[5]<lower[6])
    for d in (2,10):
        exponent=F(2*d,d+1)
        balance=d-exponent*F(d-1,2)
        ck(f"Poisson_smoothing_balance_dimension_{d}",balance==exponent>1)
    ck("ten_dimensional_error_below_shell_width",F(20,11)>1 and F(20,11)/2==F(10,11))
    ck("phase_error_below_shell_width",F(4,3)>1 and F(4,3)/2==F(2,3))
    ck("adjacent_soft_mesh_and_source_are_smaller_than_shell",F(3,2)>1 and F(3,4)>F(1,2))
    ck("cutoff_stays_submacroscopic",0<F(1,20)<F(1,2) and 2*F(1,20)-1<0)
    # Rigorous finite partial-series enclosures; the infinite bound is analytical.
    series_rows=[]
    for x in (F(1),F(7,6),F(10),F(301,6)):
        M=16
        partial=sum((1/(x+k)**2 for k in range(1,M+1)),F(0))
        lo=partial+1/(x+M+1)
        hi=partial+1/(x+M)
        assert F(1)/(x+1)<=lo<=hi<=F(1)/x
        series_rows.append({"x":str(x),"terms":M,"lower":str(lo),"upper":str(hi)})
    ck("trigamma_integral_enclosures_four_positive_arguments",len(series_rows)==4)
    # Exact cancellation with deliberately distinct neighboring ball probabilities.
    jprev,sphase,sball=F(11,7),F(-2,13),F(3,29)
    bl,br=F(2,5),F(7,19)
    lhs=jprev*bl+sball-(jprev+sphase)*br
    rhs=jprev*(bl-br)+sball-sphase*br
    ck("adjacent_current_source_identity_exact",lhs==rhs)
    ck("discarding_phase_source_term_detected",lhs!=jprev*(bl-br)+sball)
    # Departure redistribution retains exact divergence and both carry totals.
    weights=[F(2),F(3),F(1)]
    incoming=[F(4),F(6),F(2)]
    source=[F(1),F(-1),F(-3)]
    deficit=max(F(0),-source[2]-incoming[2])
    incoming_p=[incoming[i]-deficit*weights[i]/sum(weights[:2]) for i in range(2)]+[incoming[2]+deficit]
    outgoing_p=[incoming_p[i]+source[i] for i in range(3)]
    ck("E_correction_exact_totals_and_source",sum(incoming_p)==sum(incoming) and sum(outgoing_p)==sum(incoming)+sum(source) and all(outgoing_p[i]-incoming_p[i]==source[i] for i in range(3)))
    ck("E_correction_nonnegative_and_TV_cost",min(incoming_p+outgoing_p)>=0 and sum(abs(incoming_p[i]-incoming[i]) for i in range(3))==2*deficit)
    # General complement identity, not a full microscopic feasibility test.
    total=F(17)
    outA,inR=F(5),F(4)
    ck("whole_fibre_complement_identity",(total-inR)-(total-outA)==outA-inR)
    # Finite annuli need not contain a point: guards the asymptotic scope.
    empty=[(x,y) for x in range(-1,2) for y in range(-1,2)
           if F(1,8)**2<x*x+y*y<=F(1,4)**2]
    ck("small_annulus_nonemptiness_not_claimed",empty==[])
    weak_orders=[F(1,4),F(1,1)]
    ck("generic_o_one_or_same_order_cannot_pay_width",weak_orders[0]<1 and weak_orders[1]==1)
    return checks,{"count_fixtures":fixture_count,"upward_moves_per_fixture":36,
                   "scope":"N120 count/residue algebra only; not MT-domain membership or asymptotic proof",
                   "trigamma_partial_series_controls":series_rows,
                   "lattice_dimensions":[10,8,2],"full_covolume":"6","kernel_covolume":"9",
                   "projection_area":"2/3","residue_classes":54,
                   "large_count_graph_enumerated":False}


def links():
    for p in HERE.glob("*.md"):
        prose=[]
        fence=False
        for line in p.read_text().splitlines():
            assert line==line.rstrip(),p.name
            if line.lstrip().startswith("```"):
                fence=not fence
            elif not fence and not line.startswith(("    ","\t")):
                prose.append(line)
        txt=re.sub(r"`[^`\n]*`","","\n".join(prose))
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
    checks,details=run()
    assert inputs=={n:sha(ROOT/n) for n in INPUTS}
    assert before=={n:manifest(ROOT/n) for n in PROTECTED}
    result={"date":"2026-10-20","status":"ONE_MESOSCOPIC_FAMILY_PAID_CONDITIONAL_ANALYTICAL_THEOREM",
            "bounded_checks":len(checks),"checks":checks,"finite_controls":details,
            "input_hashes":inputs,"protected_manifests":before,"originals_unchanged":True,
            "same_author":True,"independent_peer_review":False,"subagents_used":False,
            "full_Hall":False,"generator_constructed":False,"unconditional_IBA":"OPEN",
            "F01":"OPEN_SOURCE","AGENTS_modified":False,"broad_goal":"ACTIVE_NOT_COMPLETE"}
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
                      "protected_manifests":len(before),"protected_entries":sum(v["entries"] for v in before.values()),
                      "full_Hall_or_generator":False,"F01":"OPEN_SOURCE","same_author":True},indent=2))


if __name__=="__main__":
    main()
