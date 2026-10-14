#!/usr/bin/env python3
import hashlib,sys,json
from pathlib import Path
root=Path(__file__).resolve().parent
lines=(root/'MANIFEST.sha256').read_text().splitlines();bad=[];seen=set()
for line in lines:
 sha,name=line.split('  ',1);p=root/name
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=sha:bad.append(name)
 seen.add(name)
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256' and '__pycache__' not in p.parts}
missing=actual-seen
print(json.dumps({'matching':len(lines)-len(bad),'entries':len(lines),'bad':bad,'unmanifested':sorted(missing)},indent=2))
sys.exit(0 if not bad and not missing else 1)
