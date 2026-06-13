#!/usr/bin/env python3
from __future__ import annotations
import hashlib, itertools, json, math, re, subprocess
from pathlib import Path
from typing import Any
from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
GEN=ROOT/'generated'
OUT=GEN/'p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.json'
MD=GEN/'p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.md'
P2657=GEN/'p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.json'
P2658=GEN/'p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.json'
STRICT_EQUATION_SHEET=ROOT/'STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md'
STRICT_LAGRANGIAN_DRAFT=ROOT/'STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md'
BASE={'n0_nadsoliton_origin':(0,0),'n1_light_axis':(1,0),'n2_matter_axis':(0,1),'n3_observer_downstream':(1,1),'n4_compression_probe':(2,1)}
SCALES=[0.25,0.5,1.0,2.0,3.0]; LAMBDAS=[0.0,0.05,0.2,1.0]; EXPONENTS=[1.0,1.8,2.0]; QUANTA=[1,2,3]; TOL=1e-11
NEG=['anomaly_clock_source_exported','intrinsic_anomaly_coefficient_derived','uv_unit_selected_by_anomaly','typed_metric_uv_source_theorem_exported','target_independent_beta_source_exported','canonical_unit_exported','declared_anomaly_anchor_promoted_to_source','bridge_completion_exported','role_transfer_revalidated','q_w_2191_discharged','role_bearing_ltotal_reenabled','toe_closure_claimed']
def sha(path:Path): return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
def shaj(o:Any): return hashlib.sha256(json.dumps(o,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def load(path:Path): return json.loads(path.read_text()) if path.exists() else {}
def rg_count(pat:str):
 r=subprocess.run(['rg','-n',pat,'.','-g','*.py','-g','*.md','-g','*.tex','-g','*.json','-g','!fundamental_action_reconstruction/generated/**','-g','!.git/**'],cwd=REPO,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 lines=sorted(x for x in r.stdout.splitlines() if x); return {'count':len(lines),'samples':lines[:60]}
def semantic_rg_audit():
 pats={'nonhomogeneous_anomaly_content':'nonhomogeneous|anomaly|boundary.*term|discrete phase|clock-source','scale_breaking_content':'scale orbit|scale-clock|UV unit|normalization anchor|fixed clock','homogeneous_no_go_content':'homogeneous|local.*action|integer phase|tau -> a','nonclosure_guard_content':'role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|beta source'}
 return {'tool':'rg','mode':'content-first semantic audit for nonhomogeneous/anomaly clock-source candidates, not packet-name search','patterns':{k:rg_count(v) for k,v in pats.items()}}
def latest_commit_audit():
 p=subprocess.run(['git','log','-n','12','--oneline','--name-only'],cwd=REPO,text=True,stdout=subprocess.PIPE,check=True); rows=[]; cur=None
 for line in p.stdout.splitlines():
  if not line.strip(): continue
  if re.match(r'^[0-9a-f]{7,12} ',line):
   if cur: rows.append(cur)
   a,b=line.split(' ',1); cur={'sha':a,'subject':b,'files':[]}
  elif cur: cur['files'].append(line)
 if cur: rows.append(cur)
 return rows
def dist(a,b):
 ax,ay=BASE[a]; bx,by=BASE[b]; return math.hypot(ax-bx,ay-by)
def trace(scale,exp): return sum(2.0/((scale*dist(a,b))**exp) for a,b in itertools.combinations(BASE,2))
def anomaly_action(scale,exp,lam): return trace(scale,exp)+lam
def anomaly_candidate_audit():
 rows=[]; max_phase=0.0; selector_rows=[]
 for exp in EXPONENTS:
  ref=trace(1.0,exp)
  for lam in LAMBDAS:
   vals=[anomaly_action(a,exp,lam) for a in SCALES]
   homogeneous_covariance_errors=[abs(v*(a**exp)-ref) for a,v in zip(SCALES,vals)]
   for a,v in zip(SCALES,vals):
    for n in QUANTA:
     tau=2*math.pi*n/v; ph=tau*v; max_phase=max(max_phase,abs(ph-2*math.pi*n))
     rows.append({'exponent':exp,'declared_anomaly_lambda':lam,'scale':a,'action_value':v,'quantized_tau':tau,'phase_error':abs(ph-2*math.pi*n)})
   # declared absolute action equal to scale-one value appears to select scale one for lambda>0 or homogeneous too; external in all cases
   q=anomaly_action(1.0,exp,lam)
   passing=[a for a,v in zip(SCALES,vals) if abs(v-q)<TOL]
   selector_rows.append({'exponent':exp,'declared_anomaly_lambda':lam,'declared_action_quantum_from_scale_one':q,'passing_scales':passing,'selects_scale_one':passing==[1.0],'is_intrinsic_source':False,'max_homogeneous_covariance_error':max(homogeneous_covariance_errors)})
 return {'statement':'Audited nonhomogeneous affine candidates A(a)=Tr_p(L_a)+lambda break pure homogeneous covariance when lambda != 0, but integer phase quantization still only fixes tau*A. Apparent scale selection requires a declared lambda and/or declared absolute action quantum, so no intrinsic anomaly clock source is exported.','rows':rows,'selector_rows':selector_rows,'max_phase_error':max_phase,'all_integer_phase_conditions_satisfied':max_phase<TOL,'nonzero_lambda_breaks_pure_homogeneous_covariance':all(r['max_homogeneous_covariance_error']>1e-6 for r in selector_rows if r['declared_anomaly_lambda']>0),'declared_absolute_action_selectors_are_external':all(not r['is_intrinsic_source'] for r in selector_rows),'uv_unit_selected_by_intrinsic_anomaly_now':False}
def upstream_consistency():
 p57=load(P2657); p58=load(P2658)
 return {'p2657_unique_scale_not_selected':p57.get('closure_decision',{}).get('unique_scale_selected_by_integer_phase_alone') is False,'p2658_homogeneous_no_go_verified':p58.get('closure_decision',{}).get('all_homogeneous_covariances_verified') is True,'p2658_anomaly_target_left_open':'nonhomogeneous' in p58.get('closure_decision',{}).get('next_honest_step','')}
def matrix(audit):
 return [{'candidate':'affine_nonhomogeneous_trace_plus_declared_lambda','covered_by_audit':True,'uses_external_clock_or_scale_anchor':True,'passes_as_uv_unit_source_now':False,'verdict':'blocked: lambda is declared, not derived from nadsoliton dynamics'}, {'candidate':'integer_phase_with_nonhomogeneous_action','covered_by_audit':True,'uses_external_clock_or_scale_anchor':False,'passes_as_uv_unit_source_now':False,'verdict':'blocked: tau can still compensate A(a) at every audited scale'}, {'candidate':'declared_absolute_action_quantum_for_anomaly','covered_by_audit':True,'uses_external_clock_or_scale_anchor':True,'passes_as_uv_unit_source_now':False,'verdict':'blocked: selection is imported by the declared absolute quantum'}, {'candidate':'derived_boundary_anomaly_with_fixed_coefficient','covered_by_audit':False,'uses_external_clock_or_scale_anchor':False,'passes_as_uv_unit_source_now':False,'verdict':'still open as next theorem target; not exported here'}]
def closure(audit,mat):
 return {'decision':'NONHOMOGENEOUS_ANOMALY_CLOCK_SOURCE_CANDIDATE_AUDIT__NO_INTRINSIC_UV_UNIT','professorial_verdict':'P2659 tests the next admissible opening left by P2658: a finite nonhomogeneous/anomalous action candidate.  The affine anomaly breaks pure homogeneous covariance, but only because lambda is supplied as an extra absolute datum. Integer phase quantization continues to fix only tau*A, and declared absolute action selection is still an external anchor. Therefore this is a useful candidate audit, not an anomaly clock-source theorem.','next_honest_step':'Try to derive the anomaly coefficient from a nadsoliton boundary/cocycle/phase law and prove it is invariant under scale-clock compensation; otherwise promote the result to a no-go theorem for declared nonhomogeneous anchors and keep beta=1 as gauge-fixed only.','passing_uv_unit_source_candidates':[r['candidate'] for r in mat if r['passes_as_uv_unit_source_now']],'all_integer_phase_conditions_satisfied':audit['all_integer_phase_conditions_satisfied'],'nonzero_lambda_breaks_pure_homogeneous_covariance':audit['nonzero_lambda_breaks_pure_homogeneous_covariance'],'uv_unit_selected_now':False,'beta_source_exported_now':False,'role_bearing_ltotal_now':False,'toe_closure_now':False}
def write_md(payload):
 a=payload['anomaly_candidate_audit']; d=payload['closure_decision']; lines=['# P2659/S1609 nonhomogeneous anomaly clock-source candidate audit','',f"Status: `{payload['status']}`",'','## Content-first audit']
 for k,v in payload['semantic_rg_antiduplication_audit']['patterns'].items(): lines.append(f'- `{k}`: {v["count"]} hits')
 lines += ['','## Computational witness',a['statement'],f"All integer phase conditions satisfied by compensating tau? `{a['all_integer_phase_conditions_satisfied']}`.",f"Nonzero lambda breaks pure homogeneous covariance? `{a['nonzero_lambda_breaks_pure_homogeneous_covariance']}`.",f"Declared absolute action selectors are external? `{a['declared_absolute_action_selectors_are_external']}`.",f"UV unit selected by intrinsic anomaly now? `{a['uv_unit_selected_by_intrinsic_anomaly_now']}`.",'','## Source candidate matrix','| candidate | covered? | external anchor? | source now? | verdict |','| --- | ---: | ---: | ---: | --- |']
 for r in payload['source_candidate_matrix']: lines.append(f"| `{r['candidate']}` | `{r['covered_by_audit']}` | `{r['uses_external_clock_or_scale_anchor']}` | `{r['passes_as_uv_unit_source_now']}` | {r['verdict']} |")
 lines += ['','## Verdict',d['professorial_verdict'],f"Decision: `{d['decision']}`.",f"Passing UV-unit source candidates: `{d['passing_uv_unit_source_candidates']}`.",f"Beta source exported now? `{d['beta_source_exported_now']}`.",f"Role-bearing L_total now? `{d['role_bearing_ltotal_now']}`.",f"ToE closure now? `{d['toe_closure_now']}`.",'','## Next honest step',d['next_honest_step'],'','## Negative exports']
 for k,v in payload['negative_export_flags'].items(): lines.append(f'- `{k}`: `{v}`')
 MD.write_text('\n'.join(lines)+'\n')
def main():
 GEN.mkdir(exist_ok=True); aud=anomaly_candidate_audit(); mat=matrix(aud); dec=closure(aud,mat)
 payload={'status':'P2659_NONHOMOGENEOUS_ANOMALY_CLOCK_SOURCE_CANDIDATE_AUDIT_NO_FALSE_PASS','latest_commit_audit':latest_commit_audit(),'semantic_rg_antiduplication_audit':semantic_rg_audit(),'source_hashes':{'P2657':sha(P2657),'P2658':sha(P2658),'STRICT_EQUATION_SHEET':sha(STRICT_EQUATION_SHEET),'STRICT_LAGRANGIAN_DRAFT':sha(STRICT_LAGRANGIAN_DRAFT)},'upstream_consistency':upstream_consistency(),'typed_anomaly_model':{'nodes':BASE,'ontology_order':'nadsoliton -> light -> matter -> emergent observer','candidate_form':'A(a)=Tr_p(L_a)+lambda_anomaly','no_sub_nadsoliton_information_layer':True},'anomaly_candidate_audit':aud,'source_candidate_matrix':mat,'closure_decision':dec,'negative_export_flags':{k:False for k in NEG}}
 payload['payload_sha256']=shaj({k:v for k,v in payload.items() if k!='payload_sha256'}); OUT.write_text(json.dumps(payload,indent=2,sort_keys=True,ensure_ascii=False)+'\n'); write_md(payload)
 append_once(STRICT_EQUATION_SHEET,'P2659/S1609 nonhomogeneous anomaly clock-source guard','## P2659/S1609 nonhomogeneous anomaly clock-source guard\n\n`P2659/S1609` audits the next nonhomogeneous/anomalous opening left by P2658 using `A(a)=Tr_p(L_a)+lambda_anomaly`.  Nonzero declared `lambda_anomaly` breaks pure homogeneous covariance, but integer phase quantization still fixes only `tau*A(a)`, and any absolute-action selection imports the anomaly coefficient/action quantum rather than deriving it from nadsoliton dynamics.  It exports no intrinsic UV unit, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n')
 append_once(STRICT_LAGRANGIAN_DRAFT,'P2659/S1609 nonhomogeneous anomaly Ltotal guard','## P2659/S1609 nonhomogeneous anomaly Ltotal guard\n\n`P2659/S1609` does not re-open `L_total`: a nonhomogeneous anomaly term can only become source-bearing after its coefficient is derived internally from a nadsoliton boundary/cocycle/phase law and shown not to be an imported scale/action anchor.  Until then, beta-source rerun, role-transfer rerun, and selector/source discharge remain blocked.\n')
 return payload
if __name__=='__main__':
 r=main(); print(rel(OUT)); print(rel(MD))
