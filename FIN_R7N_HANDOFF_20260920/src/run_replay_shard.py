from pathlib import Path
import argparse, json, os, shutil, subprocess, tempfile, time, re, sys
ROOT=Path('/mnt/data/r7n_repo_root')
SOURCE=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'
RUNNER=ROOT/'fin_rank7_intake_review'/'run_test_file.py'
ALL=sorted((SOURCE/'tests').glob('test_*.py'))
ap=argparse.ArgumentParser(); ap.add_argument('--indices',required=True); ap.add_argument('--out',required=True)
a=ap.parse_args(); idx=[int(x) for x in a.indices.split(',') if x]
rows=[]
for i in idx:
    p=ALL[i]
    parent=Path(tempfile.mkdtemp(prefix='r7n_shard_'))
    copy_root=parent/'fin_rank7_followup'
    shutil.copytree(SOURCE,copy_root,ignore=shutil.ignore_patterns('__pycache__','.pytest_cache'))
    env=os.environ.copy(); env['PYTHONDONTWRITEBYTECODE']='1'; env['PYTHONPATH']=os.pathsep.join([str(copy_root),str(parent),str(copy_root/'src')])
    start=time.monotonic(); cmd=[sys.executable,str(RUNNER),str(copy_root/'tests'/p.name)]
    try:
        cp=subprocess.run(cmd,cwd=copy_root,env=env,text=True,capture_output=True,timeout=180)
        stderr=cp.stderr; m=re.search(r'Ran (\d+) tests?',stderr)
        row={'index':i,'file':str(p.relative_to(SOURCE)),'returncode':cp.returncode,'seconds':time.monotonic()-start,'tests':int(m.group(1)) if m else None,'stdout':cp.stdout,'stderr':stderr}
    except subprocess.TimeoutExpired as e:
        row={'index':i,'file':str(p.relative_to(SOURCE)),'returncode':None,'timeout':True,'seconds':time.monotonic()-start,'tests':None,'stdout':str(e.stdout or ''),'stderr':str(e.stderr or '')}
    rows.append(row); shutil.rmtree(parent,ignore_errors=True); print(row['file'],row['returncode'],row['tests'],round(row['seconds'],2),flush=True)
Path(a.out).write_text(json.dumps(rows,indent=2)+'\n')
if any(r['returncode']!=0 for r in rows): raise SystemExit(1)
