#!/usr/bin/env python3
"""Deterministic manifest and ZIP writer. Authorial packaging, not scientific evaluation."""
from pathlib import Path
import hashlib, json, sys, platform, zipfile
root=Path(__file__).resolve().parent
ignore={'MANIFEST.sha256','verify_package.py' if False else '_none_'}
# Include verifier too, but exclude packaging script self hashes after the present run.
(root/'ENVIRONMENT.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'model':'FIN exact parameter formulas numerically binary64; no physical calibration','purpose':'replay environment, not proof'},indent=2)+'\n')
(root/'SOURCE_INPUTS.sha256').write_text(''.join(hashlib.sha256(x.read_bytes()).hexdigest()+'  '+str(x.relative_to(root))+'\n' for x in sorted((root/'source_context').glob('*')) if x.is_file()))
for i,t in enumerate(['KI-001','KI-002','KI-003','KI-004','KI-005']):
  allinputs=[root/'source_context'/'KI_PLAN_USER.md',root/'source_context'/'FOLD_SOURCE_BINDING.json',root/'source_context'/'theory_core.py',root/'SOURCE_INPUTS.sha256']
  if i>=1:allinputs.extend([root/'KI-001'/'THERMODYNAMIC_MATCHING.md',root/'KI-001'/'HESSIAN_AND_COEFFICIENT_MAP.json'])
  if i>=2:allinputs.extend([root/'KI-002'/'REVERSIBILITY_AND_MOBILITY.md',root/'KI-002'/'DERIVED_SOFT_CLOCKS.json'])
  if i>=3:allinputs.append(root/'KI-003'/'KINETIC_CRITICAL_LIMIT.md')
  if i>=4:allinputs.append(root/'KI-004'/'SPECIFICITY_VERDICT.json')
  txt=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(root))+'\n' for p in allinputs)
  (root/t/'INPUTS.sha256').write_text(txt)
files=sorted((p for p in root.rglob('*') if p.is_file() and not p.name.startswith('.')),key=lambda p:p.as_posix())
files=[f for f in files if f.name!='MANIFEST.sha256' and '__pycache__' not in f.parts]
(root/'MANIFEST.sha256').write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(root))+'\n' for f in files))
files.append(root/'MANIFEST.sha256')
archive=Path('/mnt/data/FIN_KI_001_KI_005_FULL_HANDOFF_20261014.zip')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in sorted(files,key=lambda p:p.as_posix()):z.write(p,arcname=root.name+'/'+str(p.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    bad=z.testzip();count=len(z.namelist())
assert bad is None
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'members':count,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'zip_integrity':True},indent=2))
