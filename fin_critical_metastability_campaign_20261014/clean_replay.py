#!/usr/bin/env python3
"""Reproducibility check in a fresh temp directory; does not change archived inputs."""
from pathlib import Path
import tempfile,shutil,subprocess,sys,hashlib,json
SRC=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='fin-mt-replay-') as tmp:
 dst=Path(tmp)/'campaign'
 shutil.copytree(SRC,dst,ignore=shutil.ignore_patterns('__pycache__','*.pyc','MANIFEST.sha256','CLEAN_REPLAY_VERIFICATION.json','NUMERICAL_DIAGNOSTICS.json','INDEPENDENT_VERIFICATION.json'))
 outputs=[]
 for cmd in (['mt_analysis.py'],['independent_verify.py']):
  p=subprocess.run([sys.executable,*cmd],cwd=dst,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
  if p.returncode!=0:raise RuntimeError(p.stdout+'\n'+p.stderr)
  outputs.append(p.stdout.strip())
 checks={}
 for name in ('NUMERICAL_DIAGNOSTICS.json','INDEPENDENT_VERIFICATION.json'):
  h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
  checks[name]={'identical_sha256':h(SRC/name)==h(dst/name),'original_sha256':h(SRC/name),'replayed_sha256':h(dst/name)}
 report={'scope':'same-local-software fresh temporary directory, no external sources or independent stochastic-proof acceptance',
         'checks':checks,'byte_identical_outputs':sum(c['identical_sha256'] for c in checks.values()),'total_outputs':len(checks), 'stdout':outputs}
 (SRC/'CLEAN_REPLAY_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
 print(f'CLEAN REPLAY {report["byte_identical_outputs"]}/{report["total_outputs"]} byte-identical scientific output files')
 if report['byte_identical_outputs']!=report['total_outputs']:raise SystemExit(3)
