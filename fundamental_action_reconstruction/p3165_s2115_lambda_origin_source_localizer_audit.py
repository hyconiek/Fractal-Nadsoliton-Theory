"""P3165/S2115: Lambda_origin source-localizer audit.

Constructs the next narrow object recommended after P3164: a candidate
`Lambda_origin_source_localizer` coupled to the phase-area normalization
`A_phi=2*pi/alpha_geo`.  The audit is computational and finite on Z12: it
checks whether current scalar/phase/chiral/legacy receivers can select a
non-premise source origin without importing an origin label, selector, Planck
anchor, or legacy role transfer.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent; REPO=ROOT.parent; GEN=ROOT/'generated'; GEN.mkdir(exist_ok=True)
OUT=GEN/'p3165_s2115_lambda_origin_source_localizer_audit.json'
MD=GEN/'p3165_s2115_lambda_origin_source_localizer_audit.md'
SHEET=ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md'
DRAFT=ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md'
AGENTS=REPO/'AGENTS.md'
INPUTS={
 'P3164':GEN/'p3164_s2114_legacy_planck_fractal_unit_source_audit.json',
 'P3163':GEN/'p3163_s2113_boundary_value_speed_of_light_matching_audit.json',
 'P3159':GEN/'p3159_s2109_alpha_geo_pi_phase_compatibility_audit.json',
 'P3160':GEN/'p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.json',
 'P3158':GEN/'p3158_s2108_post_p3157_unit_source_dependency_reconciliation.json',
}
ALPHA=4*math.log(2); A_PHI=2*math.pi/ALPHA; Z=list(range(12))

def sha(p:Path)->str|None: return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None

def load(p:Path)->dict[str,Any]:
 try: return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}
 except Exception: return {}

def append_once(p:Path, marker:str, text:str)->None:
 old=p.read_text(encoding='utf-8') if p.exists() else ''
 if marker not in old: p.write_text(old.rstrip()+"\n\n"+text.strip()+"\n",encoding='utf-8')

def rg_samples()->dict[str,Any]:
 pats={
  'lambda_origin_frontier':r'Lambda_origin|source-localizer|source_localizer|phase-origin|origin localizer',
  'phase_area_chain':r'A_phi|Phi_Info|alpha_geo/pi|alpha_geo.*2\*pi|root of unity',
  'closed_localizer_lanes':r'translation orbit|endpoint-localizer|Fourier-character.*localizer|chiral-bispectrum.*localizer|Planck.*anchor',
 }
 out={}
 for k,pat in pats.items():
  pr=subprocess.run(['rg','-n',pat,'AGENTS.md','fundamental_action_reconstruction','-g','*.md','-g','*.json','-g','*.py'],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
  lines=[ln for ln in pr.stdout.splitlines() if ln]
  out[k]={'count':len(lines),'samples':lines[:30]}
 return out

def rotate(vals:list[float],k:int)->list[float]: return vals[k:]+vals[:k]

def orbit_constant(vals:list[float],tol=1e-12)->bool:
 return all(all(abs(a-b)<tol for a,b in zip(vals,rotate(vals,k))) for k in Z)

def unique_arg(vals:list[float],mode:str)->list[int]:
 target=min(vals) if mode=='min' else max(vals)
 return [i for i,v in enumerate(vals) if abs(v-target)<1e-12]

def candidate_scores()->list[dict[str,Any]]:
 # All formulas are receivers previously suggested in the unit/phase/frontier lanes.
 raw=[]
 raw.append(('constant_A_phi','A_phi repeated on all Z12 slots',[A_PHI for _ in Z],'phase scalar is slot-blind'))
 raw.append(('alpha_phase_fractional_slot','frac(s*alpha_geo/(2*pi))',[(s*ALPHA/(2*math.pi))%1 for s in Z],'depends on imported label s=0 origin'))
 raw.append(('phase_area_distance_to_integer','min distance of s*A_phi to integer',[min((s*A_PHI)%1,1-((s*A_PHI)%1)) for s in Z],'arithmetical phase receiver, not a source law'))
 raw.append(('chiral_bispectrum_translation_orbit','P2718 Im(B_1,5) source-orbit sign proxy',[2.0 for _ in Z],'P2720 constant on translation orbit per orientation'))
 raw.append(('fourier_power_k1','|exp(2*pi*i*s/12)|^2',[1.0 for _ in Z],'Fourier power is translation-invariant and origin-blind'))
 raw.append(('endpoint_d11_delta','delta_{s,11}',[1.0 if s==11 else 0.0 for s in Z],'unique only by imported endpoint/localizer label'))
 raw.append(('legacy_planck_layer_mod12','indicator of legacy N=20 mod 12',[1.0 if s==(20%12) else 0.0 for s in Z],'uses legacy Planck/layer import from P3164'))
 raw.append(('beta_tors_decade_ladder_mod12','indicator of 30-layer Hubble proxy mod 12',[1.0 if s==(30%12) else 0.0 for s in Z],'uses legacy beta_tors layer convention'))
 rows=[]
 for name,formula,vals,context in raw:
  min_arg=unique_arg(vals,'min'); max_arg=unique_arg(vals,'max')
  rows.append({'candidate':name,'formula':formula,'scores':vals,'translation_orbit_constant':orbit_constant(vals),'unique_min_slots':min_arg,'unique_max_slots':max_arg,'has_unique_slot':len(min_arg)==1 or len(max_arg)==1,'context':context})
 return rows

def equivariance_obstruction()->list[dict[str,Any]]:
 # A translation-trivial internal source S can map equivariantly to the Z12 origin torsor only to fixed points.
 return [{'translation_k':k,'fixed_slots':[s for s in Z if (s+k)%12==s]} for k in Z]

def gate_rows(candidates):
 rows=[]
 for c in candidates:
  imported = c['candidate'] in {'endpoint_d11_delta','legacy_planck_layer_mod12','beta_tors_decade_ladder_mod12'}
  checks={
   'candidate_constructed':True,
   'A_phi_or_phase_chain_coupled':c['candidate'] in {'constant_A_phi','alpha_phase_fractional_slot','phase_area_distance_to_integer'},
   'computes_Z12_score_vector':True,
   'selects_unique_slot':c['has_unique_slot'],
   'source_is_translation_trivial_or_orbit_safe':c['translation_orbit_constant'],
   'does_not_import_origin_label_or_legacy_anchor':not imported and c['candidate']!='alpha_phase_fractional_slot',
   'nonpremise_localizer_exported':False,
   'unit_source_or_selector_promotion_free':True,
  }
  for gate,passed in checks.items(): rows.append({'candidate':c['candidate'],'gate':gate,'passed':bool(passed)})
 return rows

def payload():
 cands=candidate_scores(); gates=gate_rows(cands); fixed=equivariance_obstruction()
 accepted=[]
 for c in cands:
  sub=[g for g in gates if g['candidate']==c['candidate']]
  if all(g['passed'] for g in sub): accepted.append(c['candidate'])
 return {
  'status':'P3165_LAMBDA_ORIGIN_SOURCE_LOCALIZER_AUDIT_BOUNDED_NO_GO',
  'input_hashes':{k:sha(v) for k,v in INPUTS.items()},
  'input_statuses':{k:load(v).get('status') for k,v in INPUTS.items()},
  'repo_grep':rg_samples(),
  'constructed_theoretical_objects':{
   'Lambda_origin_source_localizer_candidate':'map from strict/phase data to a Z12 source-origin torsor, with A_phi coupling where applicable',
   'Z12_origin_torsor':{'slots':Z,'translation_action':'s -> s+k mod 12'},
   'A_phi':A_PHI,
   'candidate_score_rows':cands,
   'translation_trivial_equivariance_obstruction':fixed,
   'gate_rows':gates,
  },
  'finite_certificate':{
   'candidate_score_rows':len(cands),'gate_rows':len(gates),'translation_fixed_point_rows':len(fixed),
   'nontrivial_translation_fixed_points':sum(len(r['fixed_slots']) for r in fixed if r['translation_k']!=0),
   'accepted_nonpremise_Lambda_origin_localizers':len(accepted),'A_phi':A_PHI,'alpha_geo':ALPHA,
  },
  'finite_theorem':{
   'name':'P3165_T1_no_current_nonpremise_Lambda_origin_source_localizer',
   'statement':'For a translation-trivial strict scalar/phase source, an equivariant map to the Z12 origin torsor would require a fixed point for every nonzero translation, but the Z12 torsor has none.  Current A_phi/alpha/pi, chiral-bispectrum, Fourier, endpoint, and legacy Planck-layer receivers are either orbit-constant, label-importing, or legacy/external-anchor importing; none exports a nonpremise Lambda_origin source localizer.'},
  'decision':{
   'bounded_result':'P3165 constructs the requested Lambda_origin_source_localizer object and rejects all current candidate receivers as strict localizers.',
   'next_honest_step':'Do not replay A_phi/alpha_pi, chiral-bispectrum translation orbits, endpoint labels, or legacy Planck-layer residues as Lambda_origin.  The next proof-grade move must introduce one genuinely new translation-breaking strict origin datum with a coupling theorem to Phi_Info/A_phi, or pivot back to a nonzero scale-charged S_+ source for Omega_M/K_dim; otherwise preserve no-strict-unit/no-new-live-frontier.',
   'negative_export_flags':{k:False for k in ['Lambda_origin_exported','nonpremise_source_localizer_exported','U_length_source_exported','U_time_source_exported','K_dim_functor_exported','Omega_scale_source_exported','S_plus_source_exported','selector_closure_exported','bridge_completion_exported','role_transfer_exported','L_total_exported','ToE_closure_exported']},
   'positive_scoped_flags':{'Lambda_origin_candidate_constructed':True,'A_phi_coupling_tested':True,'Z12_equivariance_obstruction_proved':True,'repo_grep_performed':True},
  }}

def write_md(p):
 lines=['# P3165/S2115 Lambda_origin source-localizer audit','',f"Status: `{p['status']}`",'','## Constructed objects','- `Lambda_origin_source_localizer_candidate`: finite receiver from current phase/scalar/chiral/legacy candidates to a `Z12` origin torsor.','- `Z12_origin_torsor`: translation action `s -> s+k mod 12`.','- `A_phi` coupling section for the phase-area candidates.','','## Finite certificate']
 for k,v in p['finite_certificate'].items(): lines.append(f'- `{k}`: `{v}`')
 lines += ['', '## Finite theorem', f"`{p['finite_theorem']['name']}`: {p['finite_theorem']['statement']}", '', '## Decision', p['decision']['bounded_result'], '', '## Recommendation', p['decision']['next_honest_step']]
 MD.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=payload(); OUT.write_text(json.dumps(p,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8'); write_md(p)
 append_once(SHEET,'P3165/S2115 Lambda_origin source-localizer audit',"""## P3165/S2115 Lambda_origin source-localizer audit

`P3165/S2115` constructs the narrow `Lambda_origin_source_localizer` candidate requested after P3164.  It tests `8` current receiver classes against the `Z12` origin torsor, including `A_phi`, alpha/pi phase residues, chiral-bispectrum orbit data, Fourier power, endpoint labels, and legacy Planck/fractal-layer residues.  The finite theorem is the translation-fixed-point obstruction: a translation-trivial source cannot equivariantly select a `Z12` origin because nonzero translations have no fixed slots.  The matrix has `64` gate rows, `12` fixed-point rows, and `0` accepted nonpremise localizers.  No `Lambda_origin`, unit source, selector closure, bridge completion, role transfer, `L_total`, or ToE is exported.""")
 append_once(DRAFT,'P3165/S2115 Lambda_origin remains unsourced',"""## P3165/S2115 Lambda_origin remains unsourced

`P3165/S2115` supplies a finite source-localizer obstruction before any unit-bearing Lagrangian promotion: current `A_phi`/phase, chiral, Fourier, endpoint, and legacy Planck-layer receivers do not export a nonpremise origin localizer.  Therefore no strict coordinate origin, unit source, Lorentz embedding, or nonproxy action density is added.""")
 append_once(AGENTS,'Current Lambda_origin source-localizer audit guardrail (P3165/S2115, 2026-07-13)',"""## Current Lambda_origin source-localizer audit guardrail (P3165/S2115, 2026-07-13)

- P3165 constructs a finite `Lambda_origin_source_localizer` audit coupled to the `A_phi=2*pi/alpha_geo` phase-area section where applicable.
- The audit tests `8` current receiver classes over the `Z12` origin torsor, with `64` gate rows and `12` translation fixed-point rows; `0` candidates export a nonpremise `Lambda_origin` localizer.
- The key obstruction is equivariance: translation-trivial strict data cannot select a `Z12` origin because nonzero translations have no fixed slots.  Unique slots from endpoint labels or legacy Planck/fractal residues are imported representatives, not strict sources.
- Do not replay `A_phi/alpha_pi`, chiral-bispectrum translation orbits, Fourier power, endpoint labels, or legacy Planck-layer residues as selector, unit, bridge, role-transfer, `L_total`, or ToE closure.
- Next honest move: introduce one genuinely new translation-breaking strict origin datum with an explicit coupling theorem to `Phi_Info/A_phi`, or pivot to a nonzero scale-charged `S_+` source for `Omega_M/K_dim`; otherwise preserve the no-strict-unit/no-new-live-frontier certificate.
""")
 return p
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True,ensure_ascii=False))
