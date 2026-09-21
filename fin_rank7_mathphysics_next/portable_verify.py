#!/usr/bin/env python3
from pathlib import Path
import hashlib, json, sys
ROOT=Path(__file__).resolve().parent
man=ROOT/'MANIFEST.sha256'
if not man.exists():
    print('FAIL: MANIFEST.sha256 missing'); sys.exit(2)
bad=[]; count=0
for line in man.read_text().splitlines():
    if not line.strip(): continue
    h, rel=line.split('  ',1); p=ROOT/rel
    count+=1
    if not p.exists(): bad.append((rel,'missing')); continue
    got=hashlib.sha256(p.read_bytes()).hexdigest()
    if got!=h: bad.append((rel,'hash mismatch'))
ver=json.load(open(ROOT/'verification.json'))
tasks=json.load(open(ROOT/'TASKS.json'))
claims=json.load(open(ROOT/'CLAIM_REGISTER.json'))
print(f'manifest_files={count}')
print(f'manifest_bad={len(bad)}')
print(f'tasks={len(tasks["tasks"])} all_terminal={ver["all_tasks_terminal"]}')
print(f'claims={len(claims["claims"])} clean_replay_match={ver["clean_replay"]["all_scientific_fields_match"]}')
if bad:
    for x in bad[:20]: print('BAD',*x)
    sys.exit(1)
if len(tasks['tasks'])!=48 or not ver['all_tasks_terminal'] or not ver['clean_replay']['all_scientific_fields_match']:
    print('FAIL: ledger gate'); sys.exit(1)
print('PASS')
