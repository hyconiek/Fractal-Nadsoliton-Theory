"""Reproducible archive packer. Does NOT embed any previous research files."""
from pathlib import Path
from hashlib import sha256
from zipfile import ZipFile,ZIP_DEFLATED,ZipInfo
ROOT=Path(__file__).parent
OUT=ROOT.parent/'FIN_RI_001_RI_005_ONLY_HANDOFF_20261015.zip'
files=sorted([x for x in ROOT.rglob('*') if x.is_file() and '__pycache__' not in x.parts and x.name!='MANIFEST.sha256'])
(ROOT/'MANIFEST.sha256').write_text(''.join(f'{sha256(f.read_bytes()).hexdigest()}  {f.relative_to(ROOT)}\n' for f in files))
files.append(ROOT/'MANIFEST.sha256')
files.sort()
with ZipFile(OUT,'w',compression=ZIP_DEFLATED,compresslevel=9) as z:
    for f in files:
        name='fin_record_information_campaign_20261015/'+str(f.relative_to(ROOT))
        info=ZipInfo(name, date_time=(2026,10,15,0,0,0))
        info.compress_type=ZIP_DEFLATED
        info.external_attr=0o644<<16
        z.writestr(info,f.read_bytes(),compress_type=ZIP_DEFLATED,compresslevel=9)
with ZipFile(OUT) as z:
    assert z.testzip() is None
    assert len(z.namelist())==len(files)
print('ZIP',OUT,'bytes',OUT.stat().st_size,'members',len(files),'sha256',sha256(OUT.read_bytes()).hexdigest())
