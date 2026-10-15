"""Verify every declared package hash, stage input provenance and JSON validity.
Does not verify asymptotic proofs. The manifest itself is deliberately excluded.
"""
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
h=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
checks=[]
def test(name,ok):
    checks.append((name,bool(ok)))
    if not ok: print('FAIL',name)
for name in ['SOURCE_INPUTS.sha256']+[f'OBS-C0{i}/INPUTS.sha256' for i in range(1,6)]:
    for line in (ROOT/name).read_text().splitlines():
        expect,p=line.split('  ',1)
        target=(ROOT/name).parent/p if name!='SOURCE_INPUTS.sha256' else ROOT/p
        test(name+':'+p,target.is_file() and h(target)==expect)
for p in ROOT.rglob('*.json'):
    try:
        j=json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
        test('JSON:'+str(p.relative_to(ROOT)),True)
    except Exception as ex:test('JSON:'+str(p.relative_to(ROOT)),False)
if (ROOT/'MANIFEST.sha256').is_file():
    lines=(ROOT/'MANIFEST.sha256').read_text().splitlines()
    covered=set()
    for line in lines:
        expect,name=line.split('  ',1); covered.add(name)
        p=ROOT/name
        test('MANIFEST:'+name,p.is_file() and h(p)==expect)
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256' and '__pycache__' not in p.parts and not p.suffix=='.pyc'}
    test('MANIFEST_exhaustive',covered==actual)
summary={'passed':sum(v for _,v in checks),'total':len(checks),'all_pass':all(v for _,v in checks), 'purpose':'file-integrity/JSON validation, not proof of mathematical statements'}
print(json.dumps(summary))
if not summary['all_pass']:sys.exit(1)
