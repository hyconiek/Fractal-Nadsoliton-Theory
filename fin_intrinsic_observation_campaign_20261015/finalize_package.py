"""Deterministic-ish manifest/archive creation; file bytes frozen before archive.
Do not use this script to re-author proof or verdict documents.
"""
from pathlib import Path
import hashlib,zipfile,sys
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'FIN_OBS_C01_C05_FULL_HANDOFF_20261015.zip'
# Don't include bytecode or transient check caches.
for p in ROOT.rglob('*.pyc'):p.unlink()
for p in sorted(ROOT.rglob('*'),reverse=True):
    if p.is_dir() and p.name=='__pycache__':p.rmdir()
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name!='MANIFEST.sha256')
content=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n' for p in files)
(ROOT/'MANIFEST.sha256').write_text(content)
files.append(ROOT/'MANIFEST.sha256')
with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in sorted(files):
        name=ROOT.name+'/'+str(p.relative_to(ROOT))
        info=zipfile.ZipInfo(name,date_time=(2026,10,15,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
        info.external_attr=0o644<<16
        z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
with zipfile.ZipFile(OUT) as z:
    assert z.testzip() is None
    assert len(z.namelist())==len(files)
    root=ROOT.name+'/'
    entries={x[len(root):] for x in z.namelist()}
    assert entries=={str(p.relative_to(ROOT)) for p in files}
print('output',OUT,'bytes',OUT.stat().st_size,'files',len(files),'manifest_entries',len(files)-1,'sha256',hashlib.sha256(OUT.read_bytes()).hexdigest())
