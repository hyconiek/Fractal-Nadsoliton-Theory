from pathlib import Path
import hashlib, sys, zipfile
root=Path(__file__).resolve().parent.parent
manifest=Path(__file__).resolve().parent/'CONSOLIDATED_MANIFEST.sha256'
lines=manifest.read_text(encoding='utf-8').splitlines()
errors=[]
for line in lines:
    digest,rel=line.split('  ',1)
    f=root/rel
    if not f.is_file():errors.append(f'missing: {rel}')
    elif hashlib.sha256(f.read_bytes()).hexdigest()!=digest:errors.append(f'checksum: {rel}')
print(f'Checked {len(lines)} file hashes; problems: {len(errors)}')
for x in errors:print(x)
sys.exit(bool(errors))
