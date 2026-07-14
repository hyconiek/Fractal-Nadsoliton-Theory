"""P3164/S2114: legacy Planck/fractal layer unit-source audit.

Constructs the next proof-grade U_LT-facing object requested by P3163:
an audit of whether the legacy kernel's Planck/fractal-layer bookkeeping solves
bezwymiarowosc (dimensionlessness) as an internal strict source of length/time
units.  Result: it supplies a useful external-anchor receiver and layer-index
calculator, but not an import-free U_length/U_time theorem.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent; REPO=ROOT.parent; GEN=ROOT/'generated'; GEN.mkdir(exist_ok=True)
OUT=GEN/'p3164_s2114_legacy_planck_fractal_unit_source_audit.json'
MD=GEN/'p3164_s2114_legacy_planck_fractal_unit_source_audit.md'
SHEET=ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md'
DRAFT=ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md'
AGENTS=REPO/'AGENTS.md'
INPUTS={
 'P3163':GEN/'p3163_s2113_boundary_value_speed_of_light_matching_audit.json',
 'P3162':GEN/'p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json',
 'P3158':GEN/'p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json',
 'K1':ROOT/'K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md',
 'S2':ROOT/'S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md',
 'LEGACY_TEX':REPO/'TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex',
 'GAPS':REPO/'RAPORT_LUK_DODATKOWYCH_FIN_V5_2026-03-04.md',
}
BETA_TORS=0.01
C=299_792_458.0
HBAR=1.054_571_817e-34
G=6.674_30e-11
L_PLANCK=math.sqrt(HBAR*G/C**3)
T_PLANCK=L_PLANCK/C
M_PLANCK=math.sqrt(HBAR*C/G)
SCALES=[('Planck length',1.0,L_PLANCK),('proton proxy length',1e20,L_PLANCK*1e20),('legacy proton stated 1e-15m',1e20,1e-15),('cosmic horizon proxy',1e61,L_PLANCK*1e61),('legacy Hubble stated 1e26m',1e61,1e26)]

def sha(p:Path)->str|None: return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None

def load(p:Path)->dict[str,Any]:
 try: return json.loads(p.read_text()) if p.exists() else {}
 except Exception: return {}

def append_once(p:Path, marker:str, text:str)->None:
 old=p.read_text(encoding='utf-8') if p.exists() else ''
 if marker not in old: p.write_text(old.rstrip()+"\n\n"+text.strip()+"\n",encoding='utf-8')

def rg_samples()->dict[str,Any]:
 pats={
  'legacy_planck_layers':r'Planck.*fractal layers|fractal layers.*Planck|beta\^N|betators\^N|Planck Hierarchy',
  'dimensionless_external_anchor':r'dimensionless|external dimensionless|external anchor|without external|Planck non-manual|l_P|t_P|m_P',
  'current_unit_frontier':r'U_length|U_time|U_LT|Omega_M|Omega_scale|S_\+|K_dim|speed of light',
 }
 out={}
 for k,pat in pats.items():
  pr=subprocess.run(['rg','-n',pat,'TOE_FINAL_DOCUMENTATION_RELEASE_4_4_LEGACY_FULL.tex','RAPORT_LUK_DODATKOWYCH_FIN_V5_2026-03-04.md','LISTA_LUK_DO_UZUPELNIENIA_FIN_V5.md','fundamental_action_reconstruction','AGENTS.md','-g','*.md','-g','*.tex','-g','*.json','-g','*.py'],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
  lines=[x for x in pr.stdout.splitlines() if x]
  out[k]={'count':len(lines),'samples':lines[:25]}
 return out

def layer_rows():
 rows=[]
 for name,ratio,meters in SCALES:
  N_length=math.log(ratio,1/BETA_TORS) if ratio>0 else None
  N_force=2*N_length if N_length is not None else None
  rows.append({'scale':name,'length_m':meters,'length_ratio_to_lP':ratio,'N_length_decade_layers':N_length,'N_force_squared_layers':N_force})
 return rows

def planck_import_rows():
 return [
  {'unit':'U_length=l_P','formula':'sqrt(hbar*G/c^3)','value_SI':L_PLANCK,'requires':['hbar','G','c'],'strict_internal_source_exported':False},
  {'unit':'U_time=t_P','formula':'sqrt(hbar*G/c^5)=l_P/c','value_SI':T_PLANCK,'requires':['hbar','G','c'],'strict_internal_source_exported':False},
  {'unit':'U_mass=m_P','formula':'sqrt(hbar*c/G)','value_SI':M_PLANCK,'requires':['hbar','G','c'],'strict_internal_source_exported':False},
  {'unit':'layer_step','formula':'length multiplier 1/beta_tors=100','value_SI':100.0,'requires':['beta_tors=0.01 legacy parameter'],'strict_internal_source_exported':False},
 ]

def gate_rows():
 candidates=['legacy_Planck_anchor','legacy_beta_tors_layer_ladder','strict_G_Planck_bridge_QW2198','boundary_c_receiver_P3163']
 gates=['object_constructed','computable_finite_rows','dimensionless_ratio_bookkeeping','imports_no_external_c_hbar_G','exports_U_length','exports_U_time','exports_Lorentz_light_cone','kernel_split_role_safe']
 vals={
  'object_constructed':True,'computable_finite_rows':True,'dimensionless_ratio_bookkeeping':True,
  'imports_no_external_c_hbar_G':False,'exports_U_length':False,'exports_U_time':False,
  'exports_Lorentz_light_cone':False,'kernel_split_role_safe':True}
 return [{'candidate':c,'gate':g,'passed':vals[g]} for c in candidates for g in gates]

def payload():
 gates=gate_rows(); accepted=all(g['passed'] for g in gates)
 return {
  'status':'P3164_LEGACY_PLANCK_FRACTAL_UNIT_SOURCE_AUDIT_BOUNDED_NO_GO',
  'input_hashes':{k:sha(v) for k,v in INPUTS.items()},
  'input_statuses':{k:load(v).get('status') for k,v in INPUTS.items() if v.suffix=='.json'},
  'repo_grep':rg_samples(),
  'constructed_theoretical_objects':{
   'U_LT_legacy_receiver':'formal receiver that tries U_length=l_P and U_time=t_P from legacy Planck/fractal-layer bookkeeping',
   'PlanckUnitImportDAG':planck_import_rows(),
   'FractalLayerIndexCalculator':layer_rows(),
   'gate_rows':gates,
  },
  'finite_certificate':{
   'planck_import_rows':len(planck_import_rows()),'fractal_layer_rows':len(layer_rows()),'gate_rows':len(gates),
   'accepted_import_free_U_LT_sources':0,'beta_tors':BETA_TORS,'one_layer_length_multiplier':1/BETA_TORS,
   'computed_lP_m':L_PLANCK,'computed_tP_s':T_PLANCK,'computed_mP_kg':M_PLANCK,
  },
  'finite_theorem':{
   'name':'P3164_T1_legacy_planck_layers_are_external_anchor_receivers_not_strict_unit_sources',
   'statement':'The legacy Planck/fractal-layer mechanism resolves dimensionlessness by choosing Planck units and using beta_tors=0.01 as a logarithmic layer ladder.  The construction is computable and reproduces the intended 20/30-ish layer bookkeeping, but l_P and t_P require c,hbar,G and beta_tors remains a legacy parameter; hence it does not export an import-free U_length/U_time source or strict K_dim functor.'},
  'decision':{
   'bounded_result':'Legacy solved bezwymiarowosc operationally by external Planck anchoring plus dimensionless beta_tors layer ratios, not by an internal strict unit source.',
   'next_honest_step':'Do not reuse Planck/layer bookkeeping as closure.  The next honest proof-grade move is one narrow Lambda_origin_source_localizer audit, or a genuinely new nonzero scale-charged S_+ value coupled to Omega_M/K_dim; otherwise preserve the no-strict-unit certificate.',
   'negative_export_flags':{k:False for k in ['U_length_source_exported','U_time_source_exported','U_LT_theorem_exported','K_dim_functor_exported','Omega_scale_source_exported','S_plus_source_exported','Lorentz_metric_exported','observed_light_exported','bridge_completion_exported','role_transfer_exported','L_total_exported','ToE_closure_exported']},
   'positive_scoped_flags':{'legacy_planck_receiver_audited':True,'fractal_layer_calculator_constructed':True,'repo_grep_performed':True}
  }}

def write_md(p):
 lines=['# P3164/S2114 legacy Planck/fractal unit-source audit','',f"Status: `{p['status']}`",'','## Constructed objects','- `U_LT_legacy_receiver`: tests whether legacy `l_P/t_P` anchoring sources `U_length/U_time`.','- `PlanckUnitImportDAG`: explicit `c/hbar/G` dependency table.','- `FractalLayerIndexCalculator`: finite beta-layer rows for Planck/proton/cosmic scales.','','## Finite certificate']
 for k,v in p['finite_certificate'].items(): lines.append(f'- `{k}`: `{v}`')
 lines += ['', '## Finite theorem', f"`{p['finite_theorem']['name']}`: {p['finite_theorem']['statement']}", '', '## Decision', p['decision']['bounded_result'], '', '## Recommendation', p['decision']['next_honest_step']]
 MD.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=payload(); OUT.write_text(json.dumps(p,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8'); write_md(p)
 append_once(SHEET,'P3164/S2114 legacy Planck/fractal unit-source audit',"""## P3164/S2114 legacy Planck/fractal unit-source audit

`P3164/S2114` audits how the legacy kernel lane handled dimensionlessness: it used Planck anchoring (`l_P,t_P,m_P` from `c,hbar,G`) plus the dimensionless `beta_tors=0.01` fractal-layer ladder.  The finite audit constructs `U_LT_legacy_receiver`, a `PlanckUnitImportDAG`, and a `FractalLayerIndexCalculator`; it records `4` Planck/layer import rows, `5` layer rows, and `32` gate rows.  The result is bounded no-go for strict units: legacy Planck/layer bookkeeping is an external-anchor receiver and ratio calculator, not an import-free `U_length/U_time`, `K_dim`, `Omega_scale`, Lorentz/light-cone, bridge-completion, role-transfer, `L_total`, or ToE source theorem.""")
 append_once(DRAFT,'P3164/S2114 legacy Planck units remain external receiver',"""## P3164/S2114 legacy Planck units remain external receiver

`P3164/S2114` checks the requested legacy Planck/fractal-layer route before any speed-of-light or unit promotion.  It confirms computable layer bookkeeping, but `l_P` and `t_P` import `c,hbar,G`; no strict Lagrangian/EOM unit source or Lorentzian light-cone embedding is exported.""")
 append_once(AGENTS,'Current legacy Planck/fractal unit-source audit guardrail (P3164/S2114, 2026-07-13)',"""## Current legacy Planck/fractal unit-source audit guardrail (P3164/S2114, 2026-07-13)

- P3164 audits the legacy solution to dimensionlessness: external Planck anchoring (`l_P,t_P,m_P` from `c,hbar,G`) plus `beta_tors=0.01` fractal-layer ratios.
- The audit constructs `U_LT_legacy_receiver`, `PlanckUnitImportDAG`, and `FractalLayerIndexCalculator`; the finite matrix has `4` Planck/layer import rows, `5` layer rows, and `32` gate rows.
- This is useful operational bookkeeping, not strict closure: current artifacts still do not export import-free `U_length`, `U_time`, `K_dim`, `Omega_scale`, `S_+`, Lorentz metric/light-cone embedding, bridge completion, role transfer, `L_total`, or ToE.
- Do not transfer legacy Planck/fractal role claims onto the strict kernel without a bridge and role-transfer theorem.
- Next honest move: either supply a genuinely new nonzero scale-charged `S_+` value coupled to `Omega_M/K_dim`, or run one narrow `Lambda_origin_source_localizer` audit; otherwise preserve the no-strict-unit/no-new-live-frontier certificate.
""")
 return p
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True,ensure_ascii=False))
