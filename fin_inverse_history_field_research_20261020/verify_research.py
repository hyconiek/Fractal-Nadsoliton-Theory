#!/usr/bin/env python3
"""Exact small controls and spectral arithmetic; not proof automation."""
import argparse
import ast
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations_with_replacement, product
import json
from math import comb, factorial
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    "AGENTS.md", "research_tasks/00_EXECUTION_RULES.md",
    "fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md",
    "fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md",
    "fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md",
    "fin_intrinsic_scaling_research_20261020/GLOBAL_SCALING_THEOREM.md",
    "fin_finite_history_source_research_20261020/HANDOFF.md",
    "fin_finite_history_source_research_20261020/CONTRACT.md",
    "fin_global_activity_defect_self_stress_20261020/USER_DIRECTION.md",
    "fin_global_activity_defect_self_stress_20261020/interval_stress.py",
    "fin_rank7_followup/inputs/fin_handoff_audit/results.json",
    "fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py",
    "fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json",
    "fin_critical_metastability_review_20261015/verification.json",
    "FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md",
    "FIN_ST6897_ST8396_Hundred_Round_Entropic_Frenetic_and_FDT_Report_EN.md",
]
PROTECTED = [
    "fin_hyp_history_entropy_20261020", "fin_hyp_history_entropy_review_20261020",
    "fin_reference_transport_research_20261020", "fin_intrinsic_scaling_research_20261020",
    "fin_global_activity_defect_research_20261020", "fin_global_activity_defect_self_stress_20261020",
    "fin_critical_metastability_review_20261015", "fin_uniform_critical_review_20261014",
    "fin_intrinsic_observation_review_20261015", "fin_finite_history_source_research_20261020",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen = set()
    for line in (folder/"MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert name not in seen and not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == digest, (folder, name)
        seen.add(name)
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
              if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}
    actual |= {n for n in seen if "__pycache__" in Path(n).parts or n.endswith(".pyc")}
    assert actual == seen | {"MANIFEST.sha256"}, folder
    return {"entries":len(seen),"sha256":sha(folder/"MANIFEST.sha256")}


@dataclass(frozen=True)
class Q3:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self,"a",F(self.a))
        object.__setattr__(self,"b",F(self.b))
    @staticmethod
    def coerce(x):
        return x if isinstance(x,Q3) else Q3(x)
    def __add__(self,other):
        y = self.coerce(other)
        return Q3(self.a+y.a,self.b+y.b)
    __radd__ = __add__
    def __neg__(self):
        return Q3(-self.a,-self.b)
    def __sub__(self,other):
        return self+-self.coerce(other)
    def __mul__(self,other):
        y = self.coerce(other)
        return Q3(self.a*y.a+3*self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__ = __mul__
    def __truediv__(self,other):
        other = F(other)
        return Q3(self.a/other,self.b/other)


@dataclass(frozen=True)
class I:
    lo: F
    hi: F
    def __post_init__(self):
        object.__setattr__(self,"lo",F(self.lo))
        object.__setattr__(self,"hi",F(self.hi))
        assert self.lo <= self.hi
    @staticmethod
    def coerce(x):
        return x if isinstance(x,I) else I(F(x),F(x))
    def __add__(self,other):
        y=self.coerce(other)
        return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self):
        return I(-self.hi,-self.lo)
    def __sub__(self,other):
        return self+-self.coerce(other)
    def __mul__(self,other):
        y=self.coerce(other)
        vals=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
        return I(min(vals),max(vals))
    __rmul__=__mul__
    def __truediv__(self,other):
        y=self.coerce(other)
        assert y.lo > 0 or y.hi < 0
        return self*I(1/y.hi,1/y.lo)
    def square(self):
        return I(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
    def public(self,scale=10**12):
        return {"scale":scale,"lo":int((self.lo*scale).__floor__()),"hi":int((self.hi*scale).__ceil__())}


def exp_upper(x,n=30):
    x=F(x)
    assert x>=0 and x<n+2
    partial=sum((x**k/factorial(k) for k in range(n+1)),F(0))
    return partial+(x**(n+1)/factorial(n+1))/(1-x/F(n+2))


def unit_log(x):
    assert 1<=x<=2
    y=(x-1)/(x+1)
    n=70
    lo=2*sum((y**(2*k+1)/F(2*k+1) for k in range(n)),F(0))
    tail=2*y**(2*n+1)/(F(2*n+1)*(1-y*y))
    return I(lo,lo+tail)


@lru_cache(None)
def log_box(x):
    x=F(x)
    assert x>0
    k=0
    while x<1:
        x*=2
        k-=1
    while x>=2:
        x/=2
        k+=1
    value=unit_log(x)+k*unit_log(F(2))
    scale=1<<100
    return I(F((value.lo*scale).__floor__(),scale),F((value.hi*scale).__ceil__(),scale))


def inverse(A):
    n=len(A)
    rows=[[F(x) for x in A[i]]+[F(i==j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot]=rows[pivot],rows[j]
        d=rows[j][j]
        rows[j]=[x/d for x in rows[j]]
        for i in range(n):
            if i!=j:
                d=rows[i][j]
                rows[i]=[x-d*y for x,y in zip(rows[i],rows[j])]
    return [row[n:] for row in rows]


def run():
    checks=[]
    def ck(name,ok):
        assert name not in checks and ok,name
        checks.append(name)
    r3=Q3(0,1)
    ck("sqrt3_exact_field_relation",r3*r3==Q3(3))
    cosine=[Q3(1),r3/2,Q3(F(1,2)),Q3(),Q3(F(-1,2)),-r3/2,
            Q3(-1),-r3/2,Q3(F(-1,2)),Q3(),Q3(F(1,2)),r3/2]
    ck("cosine_mass_zero",sum(cosine)==Q3())
    ck("cosine_square_norm_six",sum(x*x for x in cosine)==Q3(6))
    modes=(3,4,5,6)
    basis={k:[cosine[(k*j)%12]/(24 if k==6 else 12) for j in range(12)] for k in modes}
    ck("native_pair_energy_zero_mean",all(sum(row)==Q3() for row in basis.values()))
    for r in (1,2):
        ck(f"native_energy_missing_cosine_{r}",all(sum(basis[k][j]*cosine[(r*j)%12] for j in range(12))==Q3() for k in modes))
    expected1={(3,4):F(1,144),(4,5):F(1,144),(5,6):F(1,144)}
    expected2={(3,5):F(1,144),(4,6):F(1,144),(5,5):F(1,288)}
    for r,expected in ((1,expected1),(2,expected2)):
        coefficients={}
        for k,l in combinations_with_replacement(modes,2):
            value=F(2-(k==l),6)*sum(basis[k][j]*basis[l][j]*cosine[(r*j)%12] for j in range(12))
            if value!=Q3():
                assert value.b==0
                coefficients[(k,l)]=value.a
        ck(f"quadratic_native_cosine_{r}_polynomial_exact",coefficients==expected)
    # Orthogonality alone would fail if no low Fourier mode were actually produced.
    parity_only=[F((-1)**j) for j in range(12)]
    ck("parity_only_does_not_fake_cosine_one",sum(parity_only[j]**2*cosine[j] for j in range(12))==Q3())
    spectral=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
    geom=json.loads((ROOT/"fin_critical_metastability_review_20261015/GEOMETRY_CERTIFICATE.json").read_text())
    ck("actual_spectral_input_hash",sha(spectral)==geom["spectral_sha256"])
    mt=json.loads((ROOT/"fin_critical_metastability_review_20261015/verification.json").read_text())
    seed="fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"
    ck("versioned_A7_convention_hash",sha(ROOT/seed)==mt["input_hashes"][seed])
    raw=json.loads(spectral.read_text())["exact"]["laplacian_intervals"]
    lam={k:I(F(raw[k][0]),F(raw[k][1])) for k in modes}
    ck("required_modes_strictly_positive",all(lam[k].lo>0 for k in modes))
    g=I.coerce(12)/lam[6]
    E0=(2*lam[3]+2*lam[4]+2*lam[5]+lam[6])/24
    C1=(lam[3]*lam[4]+lam[4]*lam[5]+lam[5]*lam[6])/144
    C2=(2*lam[3]*lam[5]+2*lam[4]*lam[6]+lam[5].square())/288
    ck("critical_gain_strictly_greater_than_five",g.lo>5)
    ck("actual_quadratic_C1_above_one_tenth",C1.lo>F(1,10))
    ck("actual_quadratic_C2_strictly_positive",C2.lo>0)
    ck("whole_energy_range_paid_by_E0",(g*E0).hi<F(10,3))
    expbound=exp_upper(F(10,3))
    ck("rational_exp_Taylor_tail_below_29",expbound<29)
    field_lower=F(25,4*12*29*10)
    ck("public_C1_history_field_lower_bound",field_lower==F(5,2784) and field_lower>F(1,600))
    p0=F(1,12000)
    ck("public_Gibbs_minimum_bound",F(1,12*29**2)>p0)
    B=1+22*(24000+23999*11)
    ck("full_interval_IS_B_exact",B==6335759)
    M=24*(B+1)/p0
    ck("full_interval_entropy_curvature_bound",M==1824698880000 and M<2*10**12)
    t=F(1,10**16)
    ck("exact_finite_perturbation_positive_interval",0<t<p0/2)
    gain=t/100-10**12*t*t
    ck("finite_entropy_advantage_certified",gain==F(99,10**20) and gain>F(1,2*10**18))
    # Exact q4/16-label positive scaling fixture, not the native A7 matrix.
    q=4
    states=list(product(range(q),repeat=2))
    a=(1,2,3,2)
    aa=sum(a)
    normal=2*aa*aa-sum(x*x for x in a)
    sep=[(y-x)%q for x,y in states]
    K=[[int(x==y or sum(u!=v for u,v in zip(x,y))==1) for y in states] for x in states]
    pi=[F(a[d]*(2*aa-a[d]),q*normal) for d in sep]
    flow=[[F(K[i][j]*a[sep[i]]*a[sep[j]],q*normal) for j in range(q*q)] for i in range(q*q)]
    P=[[flow[i][j]/pi[i] for j in range(q*q)] for i in range(q*q)]
    ck("q4_exact_labelled_mass_and_rows",sum(pi)==1 and all(sum(row)==1 for row in P))
    ck("q4_exact_stationary_reversible_flows",all(sum(flow[i][j] for i in range(q*q))==pi[j] for j in range(q*q)) and all(flow[i][j]==flow[j][i] for i in range(q*q) for j in range(q*q)))
    ck("one_unmarked_self_exact",all(K[i][i]==1 for i in range(q*q)))
    ck("two_actual_labelled_continuations_per_new_separation",all(sum(K[i][j] for j in range(q*q) if sep[j]==d)==(1 if sep[i]==d else 2) for i in range(q*q) for d in range(q)))
    s=[P[i][i] for i in range(q*q)]
    ck("inverse_field_diagonal_identity",all(s[i]==F(a[sep[i]],2*aa-a[sep[i]]) for i in range(q*q)))
    ck("symmetric_transfer_ground_relation_squared",all(flow[i][j]**2==K[i][j]*pi[i]*pi[j]*s[i]*s[j] for i in range(q*q) for j in range(q*q)))
    ck("q4_series_argument_strictly_below_one",max(F(a[d]*(2*aa-a[d]),aa*aa) for d in range(q))<1 and F(1,2)<F(aa*aa,normal)<1)
    hcoeff=[F(comb(2*n,n),n*4**n) for n in range(1,5)]
    ck("positive_H_series_coefficients",hcoeff==[F(1,2),F(3,16),F(5,48),F(35,512)])
    times=[F(1,16),F(2,16),F(3,16)]
    ck("strict_curvature_at_three_fixture_levels",all(2*x*(1-x)/(1-2*x)**3>0 for x in times))
    psmall=[F(a[d]*(2*aa-a[d]),normal) for d in range(3)]
    ssmall=[F(a[d],2*aa-a[d]) for d in range(3)]
    x1,x2,x3=[log_box(x) for x in psmall]
    f1,f2,f3=[log_box(x) for x in ssmall]
    determinant=(f3-f2)*(x2-x1)-(f2-f1)*(x3-x2)
    ck("nonaffine_inverse_field_secant_certificate",determinant.lo>0)
    size=q*q
    J=[[F(i==j)+P[i][j] for j in range(size)] for i in range(size)]
    inv=inverse(J)
    ck("I_plus_P_inverse_exact",all(sum(J[i][t]*inv[t][j] for t in range(size))==F(i==j) for i in range(size) for j in range(size)))
    ck("self_probability_inverse_bound",max(sum(abs(v) for v in row) for row in inv)<=1/(2*min(s)))
    direction=[F(-1,q) if d==0 else F(1,q) if d==2 else F(0) for d in sep]
    ratio=[direction[i]/pi[i] for i in range(size)]
    logz_derivative=[sum(inv[i][j]*ratio[j] for j in range(size)) for i in range(size)]
    ck("exact_mass_zero_entropy_direction",sum(direction)==0)
    ck("scaling_linearized_margins",all(pi[i]*sum(J[i][j]*logz_derivative[j] for j in range(size))==direction[i] for i in range(size)))
    ck("scaling_normalization_derivative_zero",sum(pi[i]*logz_derivative[i] for i in range(size))==0)
    df=[2*logz_derivative[i]-ratio[i] for i in range(size)]
    curvature=-sum(direction[i]*df[i] for i in range(size))
    curvature_cap=sum(abs(v) for v in direction)*(1/min(s)+1)*max(abs(v) for v in ratio)
    ck("exact_entropy_hessian_negative_and_capped",curvature<0 and abs(curvature)<=curvature_cap)
    gradient_flow=I.coerce(0)
    for i in range(size):
        for j in range(size):
            if K[i][j]:
                dc=flow[i][j]*(logz_derivative[i]+logz_derivative[j])
                gradient_flow=gradient_flow-dc*log_box(P[i][j])
    gradient_field=-sum((direction[i]*log_box(s[i]) for i in range(size)),I.coerce(0))
    difference=gradient_flow-gradient_field
    ck("entropy_gradient_two_exact_formulations",difference.lo<=0<=difference.hi)
    # Complete-reference control: iid law has the target as its own field tilt.
    pp=[F(15,110),F(28,110),F(39,110),F(28,110)]
    complete_flow=[[p*r for r in pp] for p in pp]
    ck("bare_macro_reference_recovers_iid_energy_tilt",all(sum(row)==pp[i] for i,row in enumerate(complete_flow)) and all(complete_flow[i][j]/pp[i]==pp[j] for i in range(q) for j in range(q)))
    P_zero=[[F(K[i][j],2*q-1) for j in range(size)] for i in range(size)]
    ck("uniform_g_zero_source_has_constant_field",all(sum(row)==1 for row in P_zero)
       and all(sum(P_zero[i][j]/size for i in range(size))==F(1,size) for j in range(size))
       and {P_zero[i][i] for i in range(size)}=={F(1,7)})
    certificate={"date":"2026-10-20","status":"EXACT_ARITHMETIC_UNDER_INHERITED_SPECTRAL_BOXES",
                 "spectral_sha256":sha(spectral),"gain":"12/Lambda6; finite N2 diagnostic, not MT domain",
                 "native_coefficient_boxes":{"g":g.public(),"C1_E_squared":C1.public(),"C2_E_squared":C2.public(),"gE0":(g*E0).public()},
                 "exp_10_over_3_upper_less_than":29,"exp_upper_exact_rational":str(expbound),
                 "public_history_field_C1_strict_lower":"1/600","intermediate_strict_lower":str(field_lower),
                 "p_min_strict_lower":str(p0),"perturbed_adjacency_ratio_upper":24000,
                 "IS_B":B,"entropy_second_derivative_absolute_upper":str(M),
                 "chosen_perturbation":str(t),"entropy_advantage_strict_lower":"1/2000000000000000000",
                 "Taylor_lower":str(gain),"units":"nats per labelled intrinsic step; dimensionless statistical model",
                 "full_A7_history_matrix_numerically_solved":False,"independent_peer_review":False}
    diagnostics={"q4_labelled_states":size,"q4_exact_pi_separation":[str(x) for x in pp],
                 "q4_secant_gap_outward_box":determinant.public(),"q4_entropy_curvature":str(curvature),
                 "q4_gradient_difference_contains_zero":True,"A7_large_domain_enumerated":False}
    return checks,certificate,diagnostics


def links():
    for p in HERE.glob("*.md"):
        lines=[]
        fence=False
        for line in p.read_text().splitlines():
            assert line==line.rstrip(),p.name
            if line.lstrip().startswith("```"):
                fence=not fence
            elif not fence and not line.startswith(("    ","\t")):
                lines.append(line)
        txt=re.sub(r"`[^`\n]*`","","\n".join(lines))
        for target in re.findall(r"\[[^\]\n]+\]\(([^)\n]+)\)",txt):
            target=target.strip("<>").split("#")[0]
            if target and not target.startswith(("http:","https:")):
                assert (p.parent/target).exists(),(p.name,target)


def main():
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--record",action="store_true")
    group.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    if args.check:
        manifest(HERE)
    inputs={n:sha(ROOT/n) for n in INPUTS}
    protected={n:manifest(ROOT/n) for n in PROTECTED}
    checks,certificate,diagnostics=run()
    assert inputs=={n:sha(ROOT/n) for n in INPUTS}
    assert protected=={n:manifest(ROOT/n) for n in PROTECTED}
    result={"date":"2026-10-20","status":"CONDITIONAL_INVERSE_FIELD_AND_EXPLICIT_ONE_MEAN_OBSTRUCTION",
            "bounded_checks":len(checks),"checks":checks,"diagnostics":diagnostics,
            "input_hashes":inputs,"protected_manifests":protected,"sources_unchanged":True,
            "same_author":True,"independent_peer_review":False,"subagents_used":False,
            "new_physical_postulate_adopted":False,"strict_F01":"OPEN_SOURCE",
            "AGENTS_modified":False,"broad_goal":"ACTIVE_NOT_COMPLETE"}
    for p in HERE.glob("*.json"):
        json.loads(p.read_text())
    for p in HERE.glob("*.py"):
        ast.parse(p.read_text())
    input_text="".join(f"{h}  {n}\n" for n,h in inputs.items())
    if args.record:
        (HERE/"CERTIFICATE.json").write_text(json.dumps(certificate,indent=2)+"\n")
        (HERE/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
        (HERE/"INPUTS.sha256").write_text(input_text)
        links()
        files=sorted(p for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.sha256" and "__pycache__" not in p.parts)
        (HERE/"MANIFEST.sha256").write_text("".join(f"{sha(p)}  {p.relative_to(HERE).as_posix()}\n" for p in files))
    else:
        assert json.loads((HERE/"RESULTS.json").read_text())==result
        assert json.loads((HERE/"CERTIFICATE.json").read_text())==certificate
        assert (HERE/"INPUTS.sha256").read_text()==input_text
        links()
    if args.output:
        assert Path("/tmp") in args.output.resolve().parents
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"bounded_checks":len(checks),
                      "inputs":len(inputs),"protected_manifests":len(protected),
                      "protected_entries":sum(v["entries"] for v in protected.values()),
                      "C1_field_strict_lower":"1/600","finite_entropy_advantage_strict_lower":"1/2000000000000000000",
                      "same_author":True,"originals_or_AGENTS_changed":False},indent=2))


if __name__=="__main__":
    main()
