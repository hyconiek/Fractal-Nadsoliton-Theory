"""P3167/S2117: S_+ monomial source exhaustion.

After P3166 closes the finite binary-origin datum class, the remaining honest
unit-source frontier is the independent S_+ / Omega_M / K_dim scale-charged
source obligation.  This audit constructs a finite monomial candidate family
from current positive dimensionless strict/legacy receivers and tests whether
any monomial can be a weight-one R_{>0}-equivariant source datum.

Result: many positive nonzero receivers exist, but every monomial of weight-zero
inputs remains weight zero.  Formal weight-one carriers such as Omega_M have the
right representation but no nonzero strict source value, so they are circular.
"""
from __future__ import annotations

import hashlib, itertools, json, math, subprocess
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent; REPO=ROOT.parent; GEN=ROOT/'generated'; GEN.mkdir(exist_ok=True)
OUT=GEN/'p3167_s2117_s_plus_monomial_source_exhaustion.json'
MD=GEN/'p3167_s2117_s_plus_monomial_source_exhaustion.md'
SHEET=ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md'
DRAFT=ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md'
AGENTS=REPO/'AGENTS.md'
INPUTS={
 'P3166':GEN/'p3166_s2116_binary_origin_datum_exhaustive_intake.json',
 'P3162':GEN/'p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json',
 'P3161':GEN/'p3161_s2111_omega_scale_positive_torsor_source_law_audit.json',
 'P3158':GEN/'p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json',
 'P3157':GEN/'p3157_s2107_omega_dim_mass_unit_torsor_audit.json',
}
ALPHA=4*math.log(2); A_PHI=2*math.pi/ALPHA; BETA_TORS=0.01; Z12_GAP=2-2*math.cos(2*math.pi/12)
BASIS=[
 ('alpha_geo',ALPHA,0,'strict dimensionless amplitude/shape normalization'),
 ('A_phi',A_PHI,0,'dimensionless phase-area section'),
 ('entropy_count_16',math.exp(ALPHA),0,'dimensionless entropy/cardinality count'),
 ('beta_tors_legacy',BETA_TORS,0,'legacy dimensionless damping ratio; not role-transferred'),
 ('strict_beta_receiver',1.0,0,'strict damping receiver normalized to one'),
 ('z12_laplacian_gap',Z12_GAP,0,'finite Z12 spectral gap, dimensionless'),
]
EXP_RANGE=range(-2,3)

def sha(p:Path)->str|None: return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None

def load(p:Path)->dict[str,Any]:
 try: return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}
 except Exception: return {}

def append_once(p:Path, marker:str, text:str)->None:
 old=p.read_text(encoding='utf-8') if p.exists() else ''
 if marker not in old: p.write_text(old.rstrip()+"\n\n"+text.strip()+"\n",encoding='utf-8')

def rg_samples():
 pats={
  's_plus_frontier':r'S_\+|Omega_M|Omega_scale|K_dim|positive torsor|scale-charged|mass-unit',
  'dimensionless_receivers':r'alpha_geo|A_phi|entropy_count|beta_tors|z12_laplacian|dimensionless',
  'closed_unit_routes':r'Planck|apparatus|selector|bridge completion|role transfer|no-strict-unit',
 }
 out={}
 for k,pat in pats.items():
  pr=subprocess.run(['rg','-n',pat,'AGENTS.md','fundamental_action_reconstruction','-g','*.md','-g','*.json','-g','*.py'],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
  lines=[ln for ln in pr.stdout.splitlines() if ln]
  out[k]={'count':len(lines),'samples':lines[:30]}
 return out

def monomial_rows():
 rows=[]; positives=0; weight_one=0
 for exps in itertools.product(EXP_RANGE, repeat=len(BASIS)):
  if all(e==0 for e in exps): continue
  val=1.0; weight=0; factors=[]
  for (name,x,w,_),e in zip(BASIS,exps):
   if e:
    val*=x**e; weight+=w*e; factors.append(f'{name}^{e}')
  positive=math.isfinite(val) and val>0
  if positive: positives+=1
  if weight==1: weight_one+=1
  rows.append({'exponents':dict((BASIS[i][0],exps[i]) for i in range(len(BASIS)) if exps[i]),'formula':' * '.join(factors),'value':val,'R_plus_weight':weight,'positive_nonzero':positive,'weight_one_candidate':weight==1,'accepted_S_plus_source':False,'blocker':'all basis inputs have R_+ weight 0, so monomial weight is 0'})
 return rows, positives, weight_one

def formal_rows():
 return [
  {'candidate':'Omega_M','R_plus_weight':1,'positive_nonzero_source_value_exported':False,'couples_to_Omega_M_K_dim':True,'accepted_S_plus_source':False,'blocker':'formal torsor representative; using it as source is circular'},
  {'candidate':'U_readout','R_plus_weight':1,'positive_nonzero_source_value_exported':False,'couples_to_Omega_M_K_dim':True,'accepted_S_plus_source':False,'blocker':'placeholder/readout unit, not strict nadsoliton source'},
  {'candidate':'sqrt_mu2_Higgs','R_plus_weight':1,'positive_nonzero_source_value_exported':False,'couples_to_Omega_M_K_dim':True,'accepted_S_plus_source':False,'blocker':'imports unsourced Higgs mass parameter'},
 ]

def gate_rows(monomial_count:int):
 candidates=['dimensionless_monomial_family','formal_weight_one_carriers','Planck_apparatus_imports']
 gates={
  'candidate_family_constructed':True,
  'finite_exhaustion_or_explicit_inventory':True,
  'positive_nonzero_values_exist':True,
  'R_plus_weight_one_exported':False,
  'noncircular_strict_source_value_exported':False,
  'couples_to_Omega_M_K_dim_without_import':False,
  'avoids_Planck_apparatus_selector_bridge_imports':True,
  'accepted_S_plus_source':False,
 }
 return [{'candidate_class':c,'gate':g,'passed':v,'monomial_rows':monomial_count if c=='dimensionless_monomial_family' else None} for c in candidates for g,v in gates.items()]

def payload():
 rows, positives, weight_one=monomial_rows(); formals=formal_rows(); gates=gate_rows(len(rows))
 sample_extremes=sorted(rows, key=lambda r:r['value'])[:5]+sorted(rows, key=lambda r:r['value'])[-5:]
 return {
  'status':'P3167_S_PLUS_MONOMIAL_SOURCE_EXHAUSTION_BOUNDED_NO_GO',
  'input_hashes':{k:sha(v) for k,v in INPUTS.items()},
  'input_statuses':{k:load(v).get('status') for k,v in INPUTS.items()},
  'repo_grep':rg_samples(),
  'constructed_theoretical_objects':{
   'S_plus_monomial_candidate_family':'integer-exponent monomials over current positive dimensionless strict/legacy receivers',
   'basis_receivers':[{'name':n,'value':v,'R_plus_weight':w,'context':c} for n,v,w,c in BASIS],
   'exponent_range':list(EXP_RANGE),
   'monomial_sample_extremes':sample_extremes,
   'formal_weight_one_carrier_rows':formals,
   'gate_rows':gates,
  },
  'finite_certificate':{
   'basis_receivers':len(BASIS),'exponent_grid_size_excluding_zero':len(rows),'positive_nonzero_monomials':positives,'weight_one_monomials':weight_one,'formal_weight_one_carriers':len(formals),'accepted_S_plus_sources':0,'gate_rows':len(gates)},
  'finite_theorem':{
   'name':'P3167_T1_dimensionless_monomials_cannot_source_S_plus',
   'statement':'Every monomial in positive dimensionless inputs has R_{>0} weight zero, so no finite product/power of alpha_geo, A_phi, entropy counts, beta_tors, strict beta, or Z12 spectral gaps can be an S_+ datum in the weight-one representation.  Formal weight-one symbols such as Omega_M have no nonzero strict source value and are circular as sources.'},
  'decision':{
   'bounded_result':'P3167 exhausts the declared monomial family and finds many positive receivers but no accepted S_+ source.',
   'next_honest_step':'Do not continue dimensionless scalar monomial, normalized ratio, Planck/apparatus, selector, or formal Omega_M self-source variants.  The next proof-grade move must either supply a genuinely scale-charged strict source value not built from weight-zero invariants, or issue a post-P3161-P3167 no-strict-unit/no-new-live-frontier certificate.',
   'negative_export_flags':{k:False for k in ['S_plus_source_exported','Omega_scale_source_exported','K_dim_functor_exported','Omega_M_fixed','unit_source_exported','Higgs_VEV_exported','EH_coupling_exported','selector_closure_exported','bridge_completion_exported','role_transfer_exported','L_total_exported','ToE_closure_exported']},
   'positive_scoped_flags':{'monomial_family_exhausted':True,'positive_dimensionless_receivers_confirmed':True,'repo_grep_performed':True},
  }}

def write_md(p):
 lines=['# P3167/S2117 S_plus monomial source exhaustion','',f"Status: `{p['status']}`",'','## Constructed objects','- `S_plus_monomial_candidate_family`: integer-exponent monomials over current positive dimensionless receivers.','- Explicit formal weight-one carrier inventory (`Omega_M`, `U_readout`, `sqrt_mu2_Higgs`).','- Gate matrix separating positive receiver existence from weight-one strict sourcehood.','','## Finite certificate']
 for k,v in p['finite_certificate'].items(): lines.append(f'- `{k}`: `{v}`')
 lines += ['', '## Finite theorem', f"`{p['finite_theorem']['name']}`: {p['finite_theorem']['statement']}", '', '## Decision', p['decision']['bounded_result'], '', '## Recommendation', p['decision']['next_honest_step']]
 MD.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=payload(); OUT.write_text(json.dumps(p,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8'); write_md(p)
 append_once(SHEET,'P3167/S2117 S_plus monomial source exhaustion',"""## P3167/S2117 S_plus monomial source exhaustion

`P3167/S2117` pivots from the closed `Lambda_origin` binary-origin class back to the independent `S_+ / Omega_M / K_dim` unit-source frontier.  It exhausts `15624` nontrivial integer-exponent monomials over `6` positive dimensionless receivers (`alpha_geo`, `A_phi`, `exp(alpha_geo)`, `beta_tors`, strict `beta`, and the `Z12` Laplacian gap).  All monomials are positive receivers but have `R_{>0}` weight `0`; `0` have weight `1`.  Formal weight-one carriers (`Omega_M`, `U_readout`, `sqrt(mu2)`) remain circular/placeholders without nonzero strict source values.  No `S_+`, `Omega_scale`, `K_dim`, unit source, Higgs VEV, EH coupling, selector closure, bridge completion, role transfer, `L_total`, or ToE is exported.""")
 append_once(DRAFT,'P3167/S2117 S_plus monomials do not add unit-bearing action',"""## P3167/S2117 S_plus monomials do not add unit-bearing action

`P3167/S2117` proves that finite monomials of current dimensionless scalar receivers remain weight-zero under the positive scale torsor.  Since no nonzero weight-one strict source value is exported, no mass/action unit, unit-bearing measure, nonproxy action density, Higgs VEV, EH coupling, or EOM term is added.""")
 append_once(AGENTS,'Current S_plus monomial source exhaustion guardrail (P3167/S2117, 2026-07-14)',"""## Current S_plus monomial source exhaustion guardrail (P3167/S2117, 2026-07-14)

- P3167 exhausts the finite monomial candidate family for `S_+`: `15624` nontrivial integer-exponent monomials over `6` positive dimensionless receivers (`alpha_geo`, `A_phi`, `exp(alpha_geo)`, `beta_tors`, strict `beta`, and the `Z12` Laplacian gap).
- Positive nonzero receivers are abundant, but every monomial has `R_{>0}` weight `0`; `0` monomials have weight `1` and `0` export an accepted `S_+` source.
- Formal weight-one carriers such as `Omega_M`, `U_readout`, and `sqrt(mu2)` have the correct representation shape but no nonzero strict source value, so using them as sources is circular/imported.
- Do not replay dimensionless scalar monomials, normalized ratios, Planck/apparatus imports, selector choices, or formal `Omega_M` self-source variants as unit, Higgs, EH, bridge, role-transfer, `L_total`, or ToE closure.
- Next honest move: supply a genuinely scale-charged strict source value not built from weight-zero invariants, or issue a post-P3161-P3167 no-strict-unit/no-new-live-frontier certificate.
""")
 return p
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True,ensure_ascii=False))
