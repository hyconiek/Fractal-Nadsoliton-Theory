#!/usr/bin/env python3
from __future__ import annotations
import csv,json,math,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import strict_records as sr
import strict_process_gate as sp
import strict_decision as sd

FIELDS=sr.FIELDS
FIX=ROOT/'fixtures';FIX.mkdir(exist_ok=True)
results=[]
def test(name,fn):
    try:
        ok,detail=fn();results.append({'name':name,'pass':bool(ok),'detail':detail})
    except Exception as e:
        results.append({'name':name,'pass':False,'detail':f'EXCEPTION {type(e).__name__}: {e}'})

def base_row(i=1,role='calibration'):
    return {
      'record_id':f'r{i}','dataset_role':role,'trajectory_id':f't{i}','replicate_id':f'p{i}',
      'timestamp_utc':f'2026-10-14T10:00:{i:02d}Z','elapsed_s':str(i/10),'device_tick':str(100+i),
      'reset_id':f'z{i}','phase_a_prepared':'0','phase_b_prepared':'1','absolute_a_read':'0','absolute_b_read':'1',
      'outcome':'1','validity_code':'VALID','gain_command':'3.0','kinetic_command':'heat_bath','control_field':'0',
      'temperature_c':'20','device_id':'dev1','readout_id':'ro1','compilation_hash':'abc123',
      'calibration_batch':'b1','latency_s':'0.001','measurement_protocol_version':'v1'
    }
def write_csv(name,rows,header=FIELDS):
    p=FIX/name
    with p.open('w',newline='',encoding='utf-8') as f:
      w=csv.DictWriter(f,fieldnames=header);w.writeheader();w.writerows(rows)
    return p

def status_fixture(name,rows,expect_pass,syntax_only=False,header=FIELDS):
    p=write_csv(name,rows,header);r=sr.audit_csv(p,require_physical_fields=not syntax_only)
    ok=(r['status'].startswith('PASS'))==expect_pass
    return ok,{'status':r['status'],'reason':r.get('reason'),'quarantine_reasons':[q['reason_codes'] for q in r.get('quarantine',[])]}

# Record-gate tests.
test('records_valid_single',lambda:status_fixture('valid.csv',[base_row()],True))
test('records_empty_bytes',lambda:(lambda r:(r['status']=='FAIL_CLOSED',r))(sr.audit_csv(FIX/'empty.bin') if (FIX/'empty.bin').exists() else (FIX/'empty.bin').write_bytes(b'') or sr.audit_csv(FIX/'empty.bin')))
test('records_header_only',lambda:status_fixture('header_only.csv',[],False))
r=base_row();r['elapsed_s']='NaN';test('records_nan_elapsed',lambda r=r:status_fixture('nan_elapsed.csv',[r],False))
r=base_row();r['latency_s']='Infinity';test('records_inf_latency',lambda r=r:status_fixture('inf_latency.csv',[r],False))
r=base_row();r['record_id']='';test('records_blank_record_id',lambda r=r:status_fixture('blank_id.csv',[r],False))
r=base_row();r['timestamp_utc']='not-a-time';test('records_bad_timestamp',lambda r=r:status_fixture('bad_timestamp.csv',[r],False))
r=base_row();r['timestamp_utc']='2026-10-14T12:00:00+02:00';test('records_nonutc_timestamp',lambda r=r:status_fixture('nonutc.csv',[r],False))
r=base_row();r['outcome']='INVALID';r['validity_code']='';test('records_invalid_missing_reason',lambda r=r:status_fixture('invalid_missing_reason.csv',[r],False))
r1=base_row(1);r2=base_row(2);r2['record_id']=r1['record_id'];test('records_duplicate_id',lambda:status_fixture('dup_id.csv',[r1,r2],False))
r1=base_row(1,'calibration');r2=base_row(2,'validation');r2['reset_id']=r1['reset_id'];test('records_reset_role_collision',lambda:status_fixture('role_collision.csv',[r1,r2],False))
r1=base_row(1,'calibration');r2=base_row(2,'validation');r2['trajectory_id']=r1['trajectory_id'];test('records_trajectory_role_collision',lambda:status_fixture('traj_collision.csv',[r1,r2],False))
r1=base_row(1,'calibration');r2=base_row(2,'validation');r2['replicate_id']=r1['replicate_id'];test('records_replicate_role_collision',lambda:status_fixture('rep_collision.csv',[r1,r2],False))
r1=base_row(1);r2=base_row(2);r2['trajectory_id']=r1['trajectory_id'];r2['reset_id']=r1['reset_id'];r2['replicate_id']=r1['replicate_id'];r2['elapsed_s']='0.01';test('records_elapsed_decrease',lambda:status_fixture('elapsed_decrease.csv',[r1,r2],False))
r1=base_row(1);r2=base_row(2);r2['trajectory_id']=r1['trajectory_id'];r2['reset_id']=r1['reset_id'];r2['replicate_id']=r1['replicate_id'];r2['timestamp_utc']='2026-10-14T09:59:59Z';test('records_timestamp_decrease',lambda:status_fixture('ts_decrease.csv',[r1,r2],False))
r1=base_row(1);r2=base_row(2);r2['trajectory_id']=r1['trajectory_id'];r2['reset_id']=r1['reset_id'];r2['replicate_id']=r1['replicate_id'];r2['device_tick']='1';test('records_tick_decrease',lambda:status_fixture('tick_decrease.csv',[r1,r2],False))
r=base_row();r['phase_a_prepared']='12';test('records_bad_q12',lambda r=r:status_fixture('bad_q12.csv',[r],False))
r=base_row();r['outcome']='2';test('records_noisy_mismatch_preserved',lambda r=r:(lambda p:(lambda x:(x['status'].startswith('PASS') and len(x['warnings'])==1,{'status':x['status'],'warnings':x['warnings']}))(sr.audit_csv(p)))(write_csv('mismatch.csv',[r])))
r=base_row();r['absolute_a_read']='MISSING';test('records_missing_physical_measurement_fails',lambda r=r:status_fixture('missing_phys.csv',[r],False))
r=base_row();r['absolute_a_read']='MISSING';test('records_missing_measurement_syntax_only_passes',lambda r=r:status_fixture('missing_syntax.csv',[r],True,syntax_only=True))
extra=FIELDS+['surprise'];r=base_row();r['surprise']='x';test('records_extra_column_rejected',lambda r=r:status_fixture('extra_col.csv',[r],False,header=extra))

# Process binding tests.
for model in sp.MODELS:
  for kin in sp.KINETICS:
    test(f'process_valid_{model}_{kin}',lambda model=model,kin=kin:(lambda r:(r['status'].startswith('PASS'),{'status':r['status'],'rowsum':r.get('generator_row_sum_max_abs')}))(sp.validate(model,kin,3.0,shortcut=(kin=='heat_bath'))))
p0=np.ones(12)/12;p0[0]=-0.1;p0[1]+=0.1
test('process_negative_p0_rejected',lambda:(lambda r:(r['status']=='FAIL_CLOSED',r['errors']))(sp.validate('FIN_A7','heat_bath',3.0,p0=p0)))
K=np.vstack([np.eye(12),np.zeros((1,12))]);K[0,0]=2
test('process_nonstochastic_K_rejected',lambda:(lambda r:(r['status']=='FAIL_CLOSED',r['errors']))(sp.validate('FIN_A7','heat_bath',3.0,K=K)))
test('process_unadmitted_model_rejected',lambda:(lambda r:(r['status']=='FAIL_CLOSED',r['errors']))(sp.validate('CUSTOM_NONCIRCULANT','heat_bath',3.0)))
test('process_reset_shortcut_wrong_kinetic_rejected',lambda:(lambda r:(r['status']=='FAIL_CLOSED',r['errors']))(sp.validate('FIN_A7','metropolis',3.0,shortcut=True)))

# Decision input validation tests.
def expect_decision_error(name,fn):
    def x():
      try:fn();return False,'accepted unexpectedly'
      except sd.DecisionInputError as e:return True,str(e)
    test(name,x)
expect_decision_error('decision_nan_n_rejected',lambda:sd.validate_counts([0]*13,'NaN'))
expect_decision_error('decision_inf_n_rejected',lambda:sd.validate_counts([0]*13,'Infinity'))
expect_decision_error('decision_fractional_n_rejected',lambda:sd.validate_counts([0]*13,'10.5'))
expect_decision_error('decision_zero_n_rejected',lambda:sd.validate_counts([0]*13,0))
expect_decision_error('decision_fractional_count_rejected',lambda:sd.validate_counts([0]*12+['1.5'],1))
expect_decision_error('decision_negative_count_rejected',lambda:sd.validate_counts([0]*12+[-1],1))
expect_decision_error('decision_sum_mismatch_rejected',lambda:sd.validate_counts([1]+[0]*12,2))
expect_decision_error('decision_alpha_zero_rejected',lambda:sd.confidence_radius_interval(10,'0'))
expect_decision_error('decision_alpha_one_rejected',lambda:sd.confidence_radius_interval(10,'1'))
expect_decision_error('decision_alpha_nan_rejected',lambda:sd.confidence_radius_interval(10,'NaN'))
# Valid threshold arithmetic.
test('decision_threshold_radius_below_rho_over_3',lambda:(lambda e:(e['upper']<=0.00119114/3, e))(sd.confidence_radius_interval(40278676,'0.025')))
test('decision_previous_n_minus_one_above_target',lambda:(lambda e:(e['lower']>0.00119114/3, e))(sd.confidence_radius_interval(40278675,'0.025')))
# Boundary-open semantics forced with a deliberately broad numerical guard.
counts=[40278676]+[0]*12
test('decision_numerical_boundary_can_fail_open',lambda:(lambda r:(r['outcome']=='NUMERICAL_BOUND_OPEN',{'outcome':r['outcome'],'boundary':r['numerical_boundary_open']}))(sd.classify(counts,40278676,numerical_guard='1')))

summary={'status':'PASS' if all(x['pass'] for x in results) else 'FAIL','passed':sum(x['pass'] for x in results),'total':len(results),'failed':[x for x in results if not x['pass']],'tests':results,
         'scope':'software/tooling checks only; no physical evidence, iid certification, or calibration coverage'}
(ROOT/'UNIT_TESTS.json').write_text(json.dumps(summary,indent=2,sort_keys=True,default=str)+'\n')
print(json.dumps({'status':summary['status'],'passed':summary['passed'],'total':summary['total'],'failed_names':[x['name'] for x in summary['failed']]},indent=2))
raise SystemExit(0 if summary['status']=='PASS' else 1)
