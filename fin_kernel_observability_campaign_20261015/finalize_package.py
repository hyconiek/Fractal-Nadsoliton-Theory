#!/usr/bin/env python3
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED,ZipInfo
import hashlib,subprocess,sys
R=Path(__file__).resolve().parent
allfiles=sorted(p for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ('MANIFEST.sha256','PACKAGE_VERIFICATION.json'))
manifest=''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(R)}\n' for p in allfiles)
(R/'MANIFEST.sha256').write_text(manifest)
subprocess.run([sys.executable,str(R/'verify_package.py')],check=True)
archive=R.parent/'FIN_KO_001_KO_005_ONLY_HANDOFF_20261015.zip'
with ZipFile(archive,'w',compression=ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(x for x in R.rglob('*') if x.is_file() and '__pycache__' not in x.parts):
  relative=R.name+'/'+str(p.relative_to(R))
  info=ZipInfo(relative,date_time=(2026,10,15,12,0,0));info.compress_type=ZIP_DEFLATED;info.external_attr=0o644 << 16
  z.writestr(info,p.read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
print('Created',archive,archive.stat().st_size,'bytes; files',len(z.namelist()))
