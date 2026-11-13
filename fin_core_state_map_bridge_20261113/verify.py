from pathlib import Path
import hashlib,json,subprocess,sys
H=Path(__file__).resolve().parent
inputs=json.loads((H/'INPUTS.json').read_text())
for p,expected in inputs['sha256'].items():
    assert hashlib.sha256((H.parent/p).read_bytes()).hexdigest()==expected,p
for p,expected in json.loads((H/'MANIFEST.json').read_text()).items():
    assert hashlib.sha256((H/p).read_bytes()).hexdigest()==expected,p
before=(H/'results.json').read_bytes()
subprocess.run([sys.executable,str(H/'check_bridge.py')],check=True,capture_output=True)
assert (H/'results.json').read_bytes()==before
r=json.loads(before);assert r['check_count']==22 and all(r['checks'].values())
print('PASS: input and output hashes; 22 exact controls; deterministic symbolic replay')
