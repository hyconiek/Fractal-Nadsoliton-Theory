#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
bad=[]; refs=0; seen=set()
for fm in sorted((ROOT/'campaigns').glob('*/FILEMAP.json')):
 data=json.loads(fm.read_text())
 for x in data['files']:
  refs+=1; p=ROOT/x['store_path']; seen.add(x['sha256'])
  if not p.exists(): bad.append(f'missing {p}')
  elif sha(p)!=x['sha256']: bad.append(f'hash {p}')
  elif p.stat().st_size!=x['size']: bad.append(f'size {p}')
for p in (ROOT/'store'/'sha256').rglob('*'):
 if p.is_file() and p.name not in seen: bad.append(f'unreferenced store object {p}')
if list(ROOT.rglob('*.zip')): bad.append('nested zip present')
# verify master manifest if present
mf=ROOT/'MASTER_MANIFEST.sha256'
if mf.exists():
 for line in mf.read_text().splitlines():
  if not line.strip(): continue
  h,rel=line.split('  ',1); p=ROOT/rel
  if not p.exists() or sha(p)!=h: bad.append(f'master manifest {rel}')
print(f'MASTER VERIFY refs={refs} unique={len(seen)} bad={len(bad)}')
for b in bad[:50]: print('BAD',b)
sys.exit(1 if bad else 0)
