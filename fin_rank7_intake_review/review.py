"""Independent repository intake: comparison, hashes, and bounded fresh-process replay.

The supplied packages are read-only inputs. Generated audit records live here.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import shutil
import tempfile
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NEW=ROOT/'fin_rank7_followup'
OLD=ROOT/'FIN_rank7_followup_handoff_bundle'

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def inventory(p):
    return {str(x.relative_to(p)):digest(x) for x in sorted(p.rglob('*'))
            if x.is_file() and '__pycache__' not in x.parts and '.pytest_cache' not in x.parts}

def compare():
    a,b=inventory(OLD/'fin_rank7_followup'),inventory(NEW)
    shared=set(a)&set(b)
    outer=[str(p.relative_to(OLD)) for p in OLD.iterdir() if p.is_file()]
    return dict(old_inner_files=len(a),new_files=len(b),
        identical=sorted(k for k in shared if a[k]==b[k]),
        changed=[dict(path=k,old_sha256=a[k],new_sha256=b[k]) for k in sorted(shared) if a[k]!=b[k]],
        absent=sorted(set(a)-set(b)),added=sorted(set(b)-set(a)),
        old_wrapper_files=sorted(outer),old_hashes=a,new_hashes=b)

def check_manifest(path,base):
    missing=[];changed=[];checked=0;invalid=[]
    for line in path.read_text().splitlines():
        if not line.strip(): continue
        m=re.fullmatch(r'([0-9a-f]{64})\s+\*?(.+)',line)
        if not m: invalid.append(line);continue
        want,name=m.groups();p=base/name
        if not p.exists(): missing.append(name)
        elif digest(p)!=want: changed.append(name)
        checked+=1
    return dict(manifest=str(path.relative_to(ROOT)),entries=checked,missing=missing,changed=changed,invalid=invalid)

def replay():
    copy_parent=Path(tempfile.mkdtemp(prefix='fin_rank7_intake_'))
    copy_root=copy_parent/'fin_rank7_followup'
    shutil.copytree(NEW,copy_root,ignore=shutil.ignore_patterns('__pycache__','.pytest_cache'))
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    env['PYTHONPATH']=os.pathsep.join([str(copy_root),str(copy_parent),str(copy_root/'src')])
    rows=[]
    for p in sorted((NEW/'tests').glob('test_*.py')):
        start=time.monotonic()
        command=[sys.executable,str(HERE/'run_test_file.py'),str(copy_root/'tests'/p.name)]
        try:
            cp=subprocess.run(command,cwd=copy_root,env=env,text=True,capture_output=True,timeout=180)
            row=dict(file=str(p.relative_to(NEW)),command=command,returncode=cp.returncode,
                     seconds=time.monotonic()-start,stdout=cp.stdout,stderr=cp.stderr)
        except subprocess.TimeoutExpired as exc:
            row=dict(file=str(p.relative_to(NEW)),command=command,returncode=None,
                     seconds=time.monotonic()-start,timeout=True,
                     stdout=str(exc.stdout or ''),stderr=str(exc.stderr or ''))
        rows.append(row)
        (HERE/'replay_unittest.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(row['file'],row['returncode'],round(row['seconds'],2),flush=True)
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    data=compare(); (HERE/'comparison.json').write_text(json.dumps(data,indent=2)+'\n')
    manifests=[check_manifest(OLD/'MANIFEST.sha256',OLD)]
    for name in ['MANIFEST.sha256','CONTINUATION_MANIFEST.sha256','FINAL_MANIFEST.sha256']:
        if (NEW/name).exists(): manifests.append(check_manifest(NEW/name,NEW))
    (HERE/'manifests.json').write_text(json.dumps(manifests,indent=2)+'\n')
    print(json.dumps(dict(old=data['old_inner_files'],new=data['new_files'],identical=len(data['identical']),
                          changed=len(data['changed']),absent=data['absent'],manifests=manifests),indent=2),flush=True)
    if args.replay: replay()

if __name__=='__main__':main()
