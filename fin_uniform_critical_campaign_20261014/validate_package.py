#!/usr/bin/env python3
"""Verify package archive and clean replay without altering producer records."""
import tempfile,shutil,subprocess,hashlib,zipfile,json,sys
from pathlib import Path
p=Path(__file__).resolve().parent
# This script is intentionally outside frozen manifest until finalization; run before manifest regen.
expected={f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ['DIAGNOSTICS.json','INDEPENDENT_VERIFICATION.json']}
with tempfile.TemporaryDirectory(prefix='fin-uf-replay-') as tmp:
 dst=Path(tmp)/p.name
 shutil.copytree(p,dst,ignore=shutil.ignore_patterns('__pycache__'))
 for cmd in ['uf_analysis.py','independent_verify.py']:
  subprocess.run([sys.executable,cmd],cwd=dst,check=True,stdout=subprocess.DEVNULL)
 found={f:hashlib.sha256((dst/f).read_bytes()).hexdigest() for f in expected}
 matches={f:expected[f]==found[f] for f in expected}
 report={'scientific_outputs_byte_equal_in_current_environment':matches,'all_equal':all(matches.values()),'replay_status':'independent fresh temp directory, same environment','proof_review':'NOT_AUTOMATED'}
 if not report['all_equal']: raise SystemExit('Numerical replay mismatch '+str(report))
 print(json.dumps(report,indent=2))
