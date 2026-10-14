#!/usr/bin/env python3
"""Fail-closed evidence CSV intake for FIN READY-001.

This module validates syntax/integrity only. It never certifies iid sampling,
physical calibration, detector validity, or provenance. Any quarantined row
blocks inferential use of the package.
"""
from __future__ import annotations
import csv, json, math, re
from collections import defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, List, Tuple

FIELDS = [
    'record_id','dataset_role','trajectory_id','replicate_id','timestamp_utc',
    'elapsed_s','device_tick','reset_id','phase_a_prepared','phase_b_prepared',
    'absolute_a_read','absolute_b_read','outcome','validity_code','gain_command',
    'kinetic_command','control_field','temperature_c','device_id','readout_id',
    'compilation_hash','calibration_batch','latency_s','measurement_protocol_version'
]
ALLOWED_ROLE = {'calibration','validation'}
ALLOWED_OUT = {str(i) for i in range(12)} | {'INVALID'}
REQUIRED_IDENTITY = {
    'record_id','dataset_role','trajectory_id','replicate_id','timestamp_utc',
    'device_tick','reset_id','device_id','readout_id','compilation_hash',
    'calibration_batch','measurement_protocol_version','kinetic_command'
}
NUMERIC_FIELDS = {'elapsed_s','latency_s','gain_command','control_field','temperature_c'}
NONNEGATIVE_FIELDS = {'elapsed_s','latency_s'}
Q12_FIELDS = {'phase_a_prepared','phase_b_prepared','absolute_a_read','absolute_b_read'}
ASCII_ID_RE = re.compile(r'^[A-Za-z0-9_.:@/+\-]+$')

class IntakeError(ValueError):
    pass


def _decimal_finite(raw: str, field: str, allow_missing: bool=False):
    if raw == '' or raw is None:
        if allow_missing:
            return None
        raise IntakeError(f'missing numeric field {field}')
    try:
        x = Decimal(raw)
    except InvalidOperation:
        raise IntakeError(f'{field}: not a decimal')
    if not x.is_finite():
        raise IntakeError(f'{field}: nonfinite decimal')
    if field in NONNEGATIVE_FIELDS and x < 0:
        raise IntakeError(f'{field}: negative')
    return x


def _utc_timestamp(raw: str):
    if not raw:
        raise IntakeError('timestamp_utc: blank')
    if raw.endswith('Z'):
        txt = raw[:-1] + '+00:00'
    else:
        txt = raw
    try:
        dt = datetime.fromisoformat(txt)
    except Exception:
        raise IntakeError('timestamp_utc: malformed ISO-8601')
    if dt.tzinfo is None:
        raise IntakeError('timestamp_utc: timezone required')
    if dt.utcoffset() != timezone.utc.utcoffset(dt):
        raise IntakeError('timestamp_utc: must be UTC')
    return dt.astimezone(timezone.utc)


def _nonblank_id(row: dict, field: str):
    v = row.get(field, '')
    if not v or not v.strip():
        raise IntakeError(f'{field}: blank')
    if len(v) > 256:
        raise IntakeError(f'{field}: too long')
    if field not in {'compilation_hash'} and not ASCII_ID_RE.match(v):
        raise IntakeError(f'{field}: unsupported characters')
    return v


def _q12_or_missing(raw: str, field: str):
    if raw in ('', 'MISSING'):
        return None
    if not raw.isdigit() or not (0 <= int(raw) <= 11):
        raise IntakeError(f'{field}: not q12 label or MISSING')
    return int(raw)


def _row_checks(row: dict, line: int, require_physical_fields: bool):
    reasons = []
    warnings = []
    parsed = {}
    try:
        for f in REQUIRED_IDENTITY:
            parsed[f] = _nonblank_id(row, f)
        if row['dataset_role'] not in ALLOWED_ROLE:
            raise IntakeError('dataset_role: not calibration/validation')
        if row['outcome'] not in ALLOWED_OUT:
            raise IntakeError('outcome: not one of 0..11/INVALID')
        if row['outcome'] == 'INVALID' and not row['validity_code'].strip():
            raise IntakeError('INVALID requires validity_code')
        parsed['timestamp_utc_dt'] = _utc_timestamp(row['timestamp_utc'])
        for f in NUMERIC_FIELDS:
            parsed[f] = _decimal_finite(row[f], f, allow_missing=(f in {'gain_command','control_field','temperature_c'}))
        tick = row['device_tick']
        if not tick.isdigit():
            raise IntakeError('device_tick: nonnegative integer required')
        parsed['device_tick_int'] = int(tick)
        for f in Q12_FIELDS:
            parsed[f] = _q12_or_missing(row[f], f)
        if require_physical_fields:
            for f in Q12_FIELDS:
                if parsed[f] is None:
                    raise IntakeError(f'{f}: required for calibration-interface evidence')
        if parsed['absolute_a_read'] is not None and parsed['absolute_b_read'] is not None and row['outcome'] != 'INVALID':
            expected = (parsed['absolute_b_read'] - parsed['absolute_a_read']) % 12
            if expected != int(row['outcome']):
                warnings.append(f'noisy label/output mismatch preserved: pair difference {expected}, output {row["outcome"]}')
    except IntakeError as e:
        reasons.append(str(e))
    return reasons, warnings, parsed


def audit_csv(path: str|Path, require_physical_fields: bool=True) -> dict:
    path = Path(path)
    if not path.exists():
        return {'status':'FAIL_CLOSED','reason':'FILE_NOT_FOUND','rows':0,'inferential_use_allowed':False}
    raw = path.read_bytes()
    if len(raw) == 0:
        return {'status':'FAIL_CLOSED','reason':'EMPTY_EVIDENCE_INPUT','rows':0,'inferential_use_allowed':False}
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        return {'status':'FAIL_CLOSED','reason':'NON_UTF8_INPUT','rows':0,'inferential_use_allowed':False}
    reader = csv.DictReader(text.splitlines())
    header = reader.fieldnames or []
    missing = sorted(set(FIELDS)-set(header))
    extra = sorted(set(header)-set(FIELDS))
    if missing or extra:
        return {'status':'FAIL_CLOSED','reason':'SCHEMA_MISMATCH','missing_columns':missing,'extra_columns':extra,'rows':0,'inferential_use_allowed':False}
    rows = list(reader)
    if not rows:
        return {'status':'FAIL_CLOSED','reason':'EMPTY_EVIDENCE_INPUT','rows':0,'inferential_use_allowed':False}

    accepted=[]; quarantined=[]; warnings=[]
    record_ids=set()
    role_by_reset={}; traj_by_reset={}; repl_by_reset={}; role_by_traj={}; role_by_repl={}
    track_obs=defaultdict(list)
    for i,row in enumerate(rows, start=2):
        reasons, row_warn, parsed = _row_checks(row, i, require_physical_fields)
        rid=row.get('record_id','')
        if rid in record_ids and rid:
            reasons.append('duplicate record_id')
        if rid:
            record_ids.add(rid)
        if not reasons:
            role=row['dataset_role']; reset=row['reset_id']; traj=row['trajectory_id']; repl=row['replicate_id']
            if reset in role_by_reset and role_by_reset[reset] != role:
                reasons.append('reset_id reused across dataset roles')
            if reset in traj_by_reset and traj_by_reset[reset] != traj:
                reasons.append('reset_id reused across trajectories')
            if reset in repl_by_reset and repl_by_reset[reset] != repl:
                reasons.append('reset_id reused across replicate_id values')
            if traj in role_by_traj and role_by_traj[traj] != role:
                reasons.append('trajectory_id reused across dataset roles')
            if repl in role_by_repl and role_by_repl[repl] != role:
                reasons.append('replicate_id reused across dataset roles')
            role_by_reset.setdefault(reset,role); traj_by_reset.setdefault(reset,traj); repl_by_reset.setdefault(reset,repl)
            role_by_traj.setdefault(traj,role); role_by_repl.setdefault(repl,role)
            track_obs[(role,traj)].append((parsed['elapsed_s'],parsed['timestamp_utc_dt'],parsed['device_tick_int'],rid,i))
        if reasons:
            quarantined.append({'line':i,'record_id':rid or None,'reason_codes':reasons,'raw_record':row})
        else:
            accepted.append(row)
            warnings.extend([{'line':i,'record_id':rid,'warning':w} for w in row_warn])

    # Monotonicity is checked only after row-level parsing succeeded.
    monotonic_fail=[]
    for key,obs in track_obs.items():
        for prev,cur in zip(obs,obs[1:]):
            if cur[0] < prev[0]: monotonic_fail.append((cur[4],cur[3],'elapsed_s decreases within trajectory'))
            if cur[1] < prev[1]: monotonic_fail.append((cur[4],cur[3],'timestamp_utc decreases within trajectory'))
            if cur[2] < prev[2]: monotonic_fail.append((cur[4],cur[3],'device_tick decreases within trajectory'))
    if monotonic_fail:
        bad_ids={r for _,r,_ in monotonic_fail}
        new_accepted=[]
        byid={x['record_id']:x for x in accepted}
        for row in accepted:
            if row['record_id'] in bad_ids:
                rs=[reason for _,rid,reason in monotonic_fail if rid==row['record_id']]
                quarantined.append({'line':next(line for line,rid,_ in monotonic_fail if rid==row['record_id']), 'record_id':row['record_id'],'reason_codes':rs,'raw_record':row})
            else:new_accepted.append(row)
        accepted=new_accepted

    package_fail = len(quarantined)>0
    roles={r:sum(x['dataset_role']==r for x in accepted) for r in sorted(ALLOWED_ROLE)}
    return {
        'status':'FAIL_CLOSED' if package_fail else 'PASS_SYNTAX_INTEGRITY_ONLY',
        'rows_total':len(rows), 'rows_accepted':len(accepted), 'rows_quarantined':len(quarantined),
        'roles_accepted':roles, 'warnings':warnings, 'quarantine':quarantined,
        'inferential_use_allowed':False if package_fail else True,
        'iid_certified':False, 'calibration_coverage_certified':False,
        'physical_provenance_certified':False,
        'note':'PASS authorizes only downstream audit, never iid/calibration/physical claims.'
    }


def main(argv=None):
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--syntax-only',action='store_true',help='allow MISSING q12 measurement fields; still no inference certificate')
    ap.add_argument('--json-out')
    ns=ap.parse_args(argv)
    result=audit_csv(ns.csv, require_physical_fields=not ns.syntax_only)
    out=json.dumps(result,indent=2,sort_keys=True,default=str)+'\n'
    if ns.json_out: Path(ns.json_out).write_text(out)
    print(out,end='')
    return 0 if result['status'].startswith('PASS') else 2

if __name__=='__main__':
    raise SystemExit(main())
