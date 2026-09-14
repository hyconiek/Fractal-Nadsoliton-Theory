#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
from src.schema import validate_task,validate_claim,validate_certificate
ROOT=Path(__file__).resolve().parent
errors=[]
nonreplayed=json.loads((ROOT/'NONREPLAYED_INPUTS.json').read_text()) if (ROOT/'NONREPLAYED_INPUTS.json').exists() else {}
# Inputs are immutable-by-policy: verify hashes.
hashes=json.loads((ROOT/'INPUT_HASHES.json').read_text())
for rel,meta in hashes.items():
    p=ROOT/rel
    if not p.exists():
        if rel in nonreplayed:
            print(f'NONREPLAYED_INPUT {rel}: {nonreplayed[rel]["reason"]}')
            continue
        errors.append(f'missing input {rel}'); continue
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h!=meta['sha256']: errors.append(f'hash mismatch {rel}')
# schema validation
try:
    for t in json.loads((ROOT/'TASKS.json').read_text()): validate_task(t)
except Exception as e: errors.append(f'TASKS: {e}')
try:
    for c in json.loads((ROOT/'CLAIMS.json').read_text()): validate_claim(c)
except Exception as e: errors.append(f'CLAIMS: {e}')
certdir=ROOT/'certificates'
for p in certdir.glob('*.json'):
    try: validate_certificate(json.loads(p.read_text()))
    except Exception as e: errors.append(f'{p.name}: {e}')
if errors:
    print('VERIFY_FAIL')
    for e in errors: print(e)
    sys.exit(1)
print('VERIFY_PASS')
print('tasks',len(json.loads((ROOT/'TASKS.json').read_text())))
print('claims',len(json.loads((ROOT/'CLAIMS.json').read_text())))
