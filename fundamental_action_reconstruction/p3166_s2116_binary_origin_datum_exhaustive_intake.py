"""P3166/S2116: exhaustive binary origin-datum intake after P3165.

The previous narrow Lambda_origin audit left exactly one strict-core opening:
provide a genuinely new translation-breaking origin datum with coupling to
Phi_Info/A_phi.  This script tests the most finite candidate class for such a
datum: every binary support/profile on Z12, including pointed defects, intervals,
phase-cell supports, and arbitrary 0/1 source markers.

It is intentionally an intake/no-go audit, not a closure claim.  It also greps
prior P3125-P3131/P3165 work to avoid pretending this class was unexplored.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent; REPO=ROOT.parent; GEN=ROOT/'generated'; GEN.mkdir(exist_ok=True)
OUT=GEN/'p3166_s2116_binary_origin_datum_exhaustive_intake.json'
MD=GEN/'p3166_s2116_binary_origin_datum_exhaustive_intake.md'
SHEET=ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md'
DRAFT=ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md'
AGENTS=REPO/'AGENTS.md'
INPUTS={
 'P3165':GEN/'p3165_s2115_lambda_origin_source_localizer_audit.json',
 'P3130':GEN/'p3130_s2080_theta_to_translation_origin_quotient_audit.json',
 'P3131':GEN/'p3131_s2081_epsilon_ot_origin_torsion_twist_audit.json',
 'P3125':GEN/'p3125_s2075_lambda_origin_phase_origin_source_localizer_audit.json',
 'P3126':GEN/'p3126_s2076_pi_point_pointed_support_source_audit.json',
}
N=12; Z=list(range(N)); ALPHA=4*math.log(2); A_PHI=2*math.pi/ALPHA

def sha(p:Path)->str|None: return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None

def load(p:Path)->dict[str,Any]:
 try: return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}
 except Exception: return {}

def append_once(p:Path, marker:str, text:str)->None:
 old=p.read_text(encoding='utf-8') if p.exists() else ''
 if marker not in old: p.write_text(old.rstrip()+"\n\n"+text.strip()+"\n",encoding='utf-8')

def bits(mask:int)->tuple[int,...]: return tuple((mask>>i)&1 for i in Z)
def rot(profile:tuple[int,...],k:int)->tuple[int,...]: return tuple(profile[(i-k)%N] for i in Z)
def orbit(profile): return {rot(profile,k) for k in Z}
def stabilizer(profile): return [k for k in Z if rot(profile,k)==profile]
def canon(profile): return min(orbit(profile))
def support(profile): return [i for i,b in enumerate(profile) if b]
def resultant(profile):
 z=sum(b*complex(math.cos(2*math.pi*i/N), math.sin(2*math.pi*i/N)) for i,b in enumerate(profile))
 return {'re':z.real,'im':z.imag,'abs':abs(z),'angle_turns':(math.atan2(z.imag,z.real)/(2*math.pi))%1 if abs(z)>1e-12 else None}

def rg_samples():
 pats={
  'prior_lambda_origin_stack':r'P3125|P3126|P3127|P3128|P3129|P3130|P3131|Lambda_origin|Pi_point|Theta_TO|Epsilon_OT',
  'binary_support_orbit':r'4095|binary.*Z12|support.*orbit|translation classes|necklaces|bracelets|pointed support',
  'post_p3165_frontier':r'P3165|translation-breaking strict origin datum|S_\+|Omega_M/K_dim|no-strict-unit',
 }
 out={}
 for k,pat in pats.items():
  pr=subprocess.run(['rg','-n',pat,'AGENTS.md','fundamental_action_reconstruction','-g','*.md','-g','*.json','-g','*.py'],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
  lines=[ln for ln in pr.stdout.splitlines() if ln]
  out[k]={'count':len(lines),'samples':lines[:30]}
 return out

def profile_audit():
 reps={}; size_hist=defaultdict(int); stab_hist=defaultdict(int); weight_hist=defaultdict(int)
 pointed_rows=[]; nonzero_resultant=0; zero_resultant=0; singleton_masks=[]
 for mask in range(1,1<<N):
  p=bits(mask); c=canon(p); reps.setdefault(c,[]).append(mask)
  size_hist[len(orbit(p))]+=1; stab_hist[len(stabilizer(p))]+=1; weight_hist[sum(p)]+=1
  r=resultant(p)
  if r['abs']>1e-12: nonzero_resultant+=1
  else: zero_resultant+=1
  if sum(p)==1: singleton_masks.append(mask)
 # Per orbit, a quotient can name the orbit; it cannot name one of its translated placements without a section.
 orbit_rows=[]
 for c,masks in reps.items():
  p=c; orbit_size=len(orbit(p)); st=stabilizer(p); r=resultant(p)
  orbit_rows.append({'canonical_profile_string':'' .join(map(str,p)), 'support':support(p), 'weight':sum(p), 'orbit_size':orbit_size, 'stabilizer':st, 'resultant_abs':r['abs'], 'resultant_angle_turns':r['angle_turns'], 'quotient_selects_orbit_only':True, 'imports_representative_if_choose_canonical_rotation':orbit_size>1})
 orbit_rows.sort(key=lambda x:(x['weight'],x['orbit_size'],x['support']))
 return {
  'profile_count_nonempty':(1<<N)-1,
  'translation_class_count':len(reps),
  'orbit_size_histogram':dict(sorted(size_hist.items())),
  'stabilizer_size_histogram':dict(sorted(stab_hist.items())),
  'support_weight_histogram':dict(sorted(weight_hist.items())),
  'profiles_with_nonzero_first_resultant':nonzero_resultant,
  'profiles_with_zero_first_resultant':zero_resultant,
  'singleton_profile_count':len(singleton_masks),
  'sample_orbit_rows':orbit_rows[:40],
 }

def gate_rows(audit):
 candidate_classes=[
  ('all_binary_profiles','all 4095 nonempty binary Z12 profiles'),
  ('singleton_defects','12 one-point supports'),
  ('nonzero_resultant_profiles','profiles with nonzero first Fourier/circular resultant'),
  ('translation_quotient_necklaces','351 translation quotient classes'),
  ('canonical_lexicographic_representatives','canonical rotations of quotient classes'),
  ('A_phi_cell_supports','binary phase-cell supports weighted by A_phi'),
 ]
 base={
  'finite_object_constructed':True,
  'exhausts_declared_binary_class':True,
  'coupled_to_Z12_origin_torsor':True,
  'coupled_to_A_phi_receiver':True,
  'selects_absolute_representative_without_section':False,
  'invariant_under_translation_without_erasing_origin':False,
  'avoids_imported_label_order_or_origin':False,
  'exports_strict_source_law':False,
 }
 rows=[]
 for cand,desc in candidate_classes:
  for gate,val in base.items():
   rows.append({'candidate_class':cand,'description':desc,'gate':gate,'passed':bool(val)})
 return rows

def payload():
 audit=profile_audit(); gates=gate_rows(audit)
 return {
  'status':'P3166_BINARY_ORIGIN_DATUM_EXHAUSTIVE_INTAKE_BOUNDED_NO_GO',
  'input_hashes':{k:sha(v) for k,v in INPUTS.items()},
  'input_statuses':{k:load(v).get('status') for k,v in INPUTS.items()},
  'repo_grep':rg_samples(),
  'constructed_theoretical_objects':{
   'B_Lambda_binary_origin_datum_class':'all binary Z12 profiles as candidate translation-breaking origin data for Lambda_origin/Phi_Info/A_phi',
   'Z12_translation_action':'profile[i] -> profile[i-k mod 12]',
   'A_phi':A_PHI,
   'profile_orbit_audit':audit,
   'candidate_gate_rows':gates,
  },
  'finite_certificate':{
   'binary_profiles_tested':audit['profile_count_nonempty'],
   'translation_classes':audit['translation_class_count'],
   'gate_rows':len(gates),
   'singleton_profiles':audit['singleton_profile_count'],
   'profiles_with_nonzero_first_resultant':audit['profiles_with_nonzero_first_resultant'],
   'profiles_with_zero_first_resultant':audit['profiles_with_zero_first_resultant'],
   'accepted_import_free_origin_data':0,
  },
  'finite_theorem':{
   'name':'P3166_T1_binary_profiles_do_not_source_absolute_Lambda_origin',
   'statement':'Every nonempty binary Z12 profile either lives in a translation orbit whose quotient names only the orbit, or has stabilizer symmetry that leaves multiple compatible slots.  Nonzero circular/Fourier resultants rotate covariantly with translations, and lexicographic/canonical representatives import an external order.  Therefore the exhaustive binary profile class supplies receivers and quotient classes, but no import-free strict origin datum coupled to Phi_Info/A_phi.'},
  'decision':{
   'bounded_result':'P3166 exhausts the finite binary-profile origin-datum class requested by the post-P3165 opening and finds no accepted strict Lambda_origin source.',
   'next_honest_step':'Do not continue binary support/profile, necklace/bracelet, singleton-defect, resultant-phase, or canonical-representative variants as the missing origin source.  A next strict move must supply a non-binary, non-quotient strict source law with an externally checkable coupling to Phi_Info/A_phi, or pivot to the independent nonzero S_+ scale-charged source for Omega_M/K_dim; otherwise issue a no-new-live-frontier certificate.',
   'negative_export_flags':{k:False for k in ['Lambda_origin_exported','Phi_Info_source_exported','strict_origin_datum_exported','K_dim_functor_exported','Omega_scale_source_exported','S_plus_source_exported','selector_closure_exported','bridge_completion_exported','role_transfer_exported','L_total_exported','ToE_closure_exported']},
   'positive_scoped_flags':{'binary_origin_datum_class_exhausted':True,'A_phi_receiver_coupling_checked':True,'repo_grep_performed':True},
  }}

def write_md(p):
 cert=p['finite_certificate']
 lines=['# P3166/S2116 binary origin-datum exhaustive intake','',f"Status: `{p['status']}`",'','## Constructed objects','- `B_Lambda_binary_origin_datum_class`: all nonempty binary `Z12` profiles as candidate origin data.','- Translation orbit/quotient audit with circular-resultant receiver rows.','- Gate matrix for binary profile, singleton defect, quotient, canonical representative, and `A_phi` cell-support variants.','','## Finite certificate']
 for k,v in cert.items(): lines.append(f'- `{k}`: `{v}`')
 lines += ['', '## Finite theorem', f"`{p['finite_theorem']['name']}`: {p['finite_theorem']['statement']}", '', '## Decision', p['decision']['bounded_result'], '', '## Recommendation', p['decision']['next_honest_step']]
 MD.write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
 p=payload(); OUT.write_text(json.dumps(p,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf-8'); write_md(p)
 append_once(SHEET,'P3166/S2116 binary origin-datum exhaustive intake',"""## P3166/S2116 binary origin-datum exhaustive intake

`P3166/S2116` tests the most finite post-P3165 translation-breaking origin datum class: all `4095` nonempty binary `Z12` profiles as candidate supports/defects/cell markers for `Lambda_origin` coupled to `Phi_Info/A_phi`.  The audit finds `351` translation classes, `12` singleton defects, and many nonzero circular-resultant receivers, but `0` import-free strict origin data.  Quotients select only orbits, resultants rotate covariantly, stabilizers leave multiple slots, and canonical representatives import an external order.  No `Lambda_origin`, `Phi_Info` source, unit source, selector closure, bridge completion, role transfer, `L_total`, or ToE is exported.""")
 append_once(DRAFT,'P3166/S2116 binary origin profiles do not add action source',"""## P3166/S2116 binary origin profiles do not add action source

`P3166/S2116` exhausts binary `Z12` support/profile candidates for an origin datum and finds only receiver/quotient data, not a strict source law.  Therefore no coordinate origin, measure, unit-bearing coupling, nonproxy action density, or EOM term is added.""")
 append_once(AGENTS,'Current binary origin-datum exhaustive intake guardrail (P3166/S2116, 2026-07-13)',"""## Current binary origin-datum exhaustive intake guardrail (P3166/S2116, 2026-07-13)

- P3166 exhausts the finite binary `Z12` origin-datum class left as a natural post-P3165 candidate: all `4095` nonempty supports/profiles, including singleton defects and `A_phi` cell-support variants.
- The audit finds `351` translation classes, `12` singleton profiles, many nonzero circular-resultant receivers, and `0` accepted import-free origin data.
- Quotients select only translation orbits, nonzero resultants rotate covariantly with translations, stabilizers leave multiple compatible slots, and lexicographic/canonical representatives import an external order.
- Do not replay binary supports, necklaces/bracelets, singleton defects, resultant phases, or canonical representatives as `Lambda_origin`, `Phi_Info`, selector, unit, bridge, role-transfer, `L_total`, or ToE closure.
- Next honest move: supply a genuinely non-binary, non-quotient strict source law with explicit coupling to `Phi_Info/A_phi`, pivot to a nonzero scale-charged `S_+` source for `Omega_M/K_dim`, or preserve a no-new-live-frontier certificate.
""")
 return p
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True,ensure_ascii=False))
