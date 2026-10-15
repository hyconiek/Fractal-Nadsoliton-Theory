from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parent
assert all((root/k).is_dir() for k in ('LI-001','LI-002','LI-003','LI-004','LI-005'))
assert json.loads((root/'INDEPENDENT_CHECKS.json').read_text())['status']=='PASS'
lines=(root/'MANIFEST.sha256').read_text().splitlines()
for ln in lines:
    sha,name=ln.split('  ',1)
    target=root/name
    assert target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest()==sha,name
print('PASS: '+str(len(lines))+' internally hashed package files; source hashes are separately recorded')
