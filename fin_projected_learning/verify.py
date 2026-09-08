"""Replay the current checkpoint; no publication generation."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import research

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

result=subprocess.run([sys.executable,'-m','unittest','discover','-s','.',
                       '-p','test_*.py','-v'],cwd=HERE,capture_output=True,text=True)
if result.returncode:
    print(result.stderr);raise SystemExit(result.returncode)
expected=json.loads((HERE/'results.json').read_text())
actual=json.loads(json.dumps(research.run()))
if actual!=expected:raise AssertionError('Checkpoint replay differs')
sources=['nadsoliton_neural_analysis.py','The FIN Kernel as an Unknown Mat.md',
         'fin_st2202_st2216_research.py','fin_replication_consistency/certify.py','AGENTS.md']
record=dict(test_count=int(re.search(r'Ran (\d+) tests',result.stderr).group(1)),
    returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
    exact_JSON_replay=True,completed_rounds=10,goal_rounds=30,
    source_sha256={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
    proof_scope='Analytic proofs plus exact spectral enclosures; numerical replay is not physical evidence.')
(HERE/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
print(f"{record['test_count']} tests passed; exact checkpoint and source replay matched; 10/30 rounds.")
