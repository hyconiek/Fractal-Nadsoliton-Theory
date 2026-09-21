from __future__ import annotations
import subprocess,time,json,re,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    files=sorted((ROOT/'tests').glob('test_*.py'))
    rows=[]; full=[]
    env=os.environ.copy(); env['PYTHONPATH']='.:..:src'
    for p in files:
        rel=str(p.relative_to(ROOT)); t=time.time()
        try:
            cp=subprocess.run([sys.executable,'-m','pytest','-q',rel],cwd=ROOT,env=env,text=True,capture_output=True,timeout=45)
            status='PASS' if cp.returncode==0 else 'FAIL'
            out=cp.stdout+cp.stderr
        except subprocess.TimeoutExpired as e:
            status='TIMEOUT_RESOURCE_STOP'; out=(e.stdout or '')+(e.stderr or '') if isinstance(e.stdout,str) else ''
            cp=None
        sec=time.time()-t
        m=re.search(r'(\d+) passed',out); passed=int(m.group(1)) if m else 0
        rows.append({'file':rel,'status':status,'passed_tests':passed,'seconds':sec})
        full.append(f'===== {rel} [{status}] {sec:.2f}s =====\n{out}\n')
        print(rel,status,passed,f'{sec:.2f}s',flush=True)
    (ROOT/'logs/P_full_sharded_tests.log').write_text('\n'.join(full))
    summary={'files':rows,'passed_files':sum(r['status']=='PASS' for r in rows),'failed_files':sum(r['status']=='FAIL' for r in rows),
             'timeout_files':sum(r['status'].startswith('TIMEOUT') for r in rows),'passed_tests':sum(r['passed_tests'] for r in rows)}
    (ROOT/'logs/P_sharded_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
