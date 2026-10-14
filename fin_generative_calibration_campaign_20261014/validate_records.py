#!/usr/bin/env python3
"""Prospective raw-record syntax audit ONLY, not proof of physical calibration/independence.
Usage: python validate_records.py input.csv
Does not impute invalid outcomes, repair timestamps, or pool trajectories as iid.
"""
import csv,sys,json
from collections import defaultdict
from pathlib import Path
FIELDS=['record_id','dataset_role','trajectory_id','replicate_id','timestamp_utc',
 'elapsed_s','device_tick','reset_id','phase_a_prepared','phase_b_prepared',
 'absolute_a_read','absolute_b_read','outcome','validity_code','gain_command',
 'kinetic_command','control_field','temperature_c','device_id','readout_id',
 'compilation_hash','calibration_batch','latency_s','measurement_protocol_version']
ALLOWED_ROLE={'calibration','validation'}
ALLOWED_OUT={str(i) for i in range(12)}|{'INVALID'}
def parse(path):
 with open(path,newline='',encoding='utf-8') as f:
  reader=csv.DictReader(f)
  missing=set(FIELDS)-set(reader.fieldnames or [])
  if missing:raise ValueError('missing columns: '+', '.join(sorted(missing)))
  rows=list(reader)
 seen=set();tracks=defaultdict(list);warnings=[]
 for idx,row in enumerate(rows,2):
  if row['record_id'] in seen:raise ValueError(f'duplicate record ID: {row["record_id"]}')
  seen.add(row['record_id'])
  if row['dataset_role'] not in ALLOWED_ROLE:raise ValueError(f'line {idx}: dataset role')
  if row['outcome'] not in ALLOWED_OUT:raise ValueError(f'line {idx}: outcome not in 13 symbols')
  if row['outcome']=='INVALID' and not row['validity_code']:raise ValueError(f'line {idx}: invalid requires reason code')
  for col in ('phase_a_prepared','phase_b_prepared','absolute_a_read','absolute_b_read'):
   if row[col] and (not row[col].isdigit() or not 0<=int(row[col])<=11):
    raise ValueError(f'line {idx}: {col} not q12 label')
  if row['absolute_a_read'] and row['absolute_b_read'] and row['outcome']!='INVALID':
   expected=(int(row['absolute_b_read'])-int(row['absolute_a_read']))%12
   if expected!=int(row['outcome']):warnings.append(f'line {idx}: physical noisy mismatch: observed output {row["outcome"]} versus pair difference {expected}; NOT overwritten')
  for col in ('elapsed_s','latency_s'):
   if row[col] and float(row[col])<0:raise ValueError(f'line {idx}: negative {col}')
  if not row['trajectory_id'] or not row['reset_id']:raise ValueError(f'line {idx}: missing trajectory/reset')
  tracks[(row['dataset_role'],row['trajectory_id'])].append((row['elapsed_s'],row['record_id']))
 # Preserved trajectories are NOT automatically independent across shots.
 return {'syntax_status':'PASS','rows':len(rows),'roles':{r:sum(x['dataset_role']==r for x in rows) for r in sorted(ALLOWED_ROLE)},
   'trajectories':len(tracks),'independent_iid_shots_certified':False,
   'warnings':warnings,'calibration_readout_clock_status':'NOT_CERTIFIED_BY_A_SYNTAX_PARSER'}
if __name__=='__main__':
 try:print(json.dumps(parse(sys.argv[1]),indent=2))
 except (ValueError,IndexError) as ex: print('INVALID:',ex,file=sys.stderr);sys.exit(2)
