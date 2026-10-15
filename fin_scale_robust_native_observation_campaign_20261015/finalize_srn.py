#!/usr/bin/env python3
"""SRN-only reproducible input hash binding, manifest validation and ZIP writer.
Never mutates or repackages prior FIN campaign inputs.
"""
from pathlib import Path
import hashlib, json, zipfile, subprocess, sys
R=Path(__file__).resolve().parent
BASE=R.parent
INPUTS=[
 BASE/'Wklejony Markdown(20261015-164447).md',
 BASE/'FIN_NOC_001_NOC_005_ONLY_HANDOFF_20261015.zip',
 BASE/'FIN_DKO_001_DKO_005_ONLY_HANDOFF_20261015.zip',
 BASE/'FIN_KO_001_KO_005_ONLY_HANDOFF_20261015.zip',
 BASE/'FIN_LI_001_LI_005_ONLY_HANDOFF_20261015.zip',
]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source=[]
for p in INPUTS:
    if p.exists():source.append(f'{sha(p)}  ../{p.name}')
(R/'SOURCE_INPUTS.sha256').write_text('\n'.join(source)+'\n')
(R/'INPUTS.sha256').write_text('\n'.join(source)+'\n')
for st in range(1,6):
    stem=f'SRN-{st:03d}'
    (R/stem/'INPUTS.sha256').write_text('\n'.join(source)+'\n')
(R/'REPLAY.md').write_text('''# SRN-only replay

From this directory run `python verify_srn.py` to regenerate `INDEPENDENT_VERIFICATION.json` and check exact rational projector identities, finite generator jumps, conditional hypergeometric fixtures, raw gauge invariance, nominal ideal angular calculations and exact-root exponent bookkeeping. Run `python finalize_srn.py` to regenerate input hashes, a comprehensive SHA-256 manifest and the SRN-only ZIP. Inputs are hashed from actual available bytes and not copied into the package.

Requirements: Python 3, numpy, sympy, mpmath. Tests are algebra/numeric diagnostics, not mechanized proofs of local Gibbs/MT or asymptotic all-N results. Read each stage's analytical proof separately. The source user plan and external archives are input files, not executed commands. Absence of a source review file must not be represented as a verified historical artifact.
''')
subprocess.run([sys.executable,str(R/'verify_srn.py')],check=True)
files=sorted(p for p in R.rglob('*') if p.is_file() and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json') and '__pycache__' not in p.parts)
(R/'MANIFEST.sha256').write_text('\n'.join(f'{sha(p)}  {p.relative_to(R)}' for p in files)+'\n')
verify=[]
for line in (R/'MANIFEST.sha256').read_text().splitlines():
    h,relative=line.split('  ',1)
    verify.append((R/relative).is_file() and sha(R/relative)==h)
meta={'manifest_entries':len(files),'manifest_matched':sum(verify),'manifest_full_coverage':len(files)==len(verify) and all(verify),'provenance_inputs_hashed':len(source),'scientific_assertions':json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())['checks_total'],'scientific_assertions_passed':json.loads((R/'INDEPENDENT_VERIFICATION.json').read_text())['checks_passed'],'scientific_scope':'bounded diagnostics, conditional analytical model theorem, not mechanically certified','archival_inputs_copied':False}
(R/'PACKAGE_VERIFICATION.json').write_text(json.dumps(meta,indent=2)+'\n')
zipdest=BASE/'FIN_SRN_001_SRN_005_ONLY_HANDOFF_20261015.zip'
with zipfile.ZipFile(zipdest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in sorted(R.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            info=zipfile.ZipInfo(str(p.relative_to(BASE)),date_time=(2026,10,15,12,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
with zipfile.ZipFile(zipdest) as z:
    z.checktest=z.testzip()
    assert z.checktest is None
    assert len(z.namelist())==len(files)+2
print('PKG',meta,'zip=',zipdest,'items=',len(files)+2,'size=',zipdest.stat().st_size)
