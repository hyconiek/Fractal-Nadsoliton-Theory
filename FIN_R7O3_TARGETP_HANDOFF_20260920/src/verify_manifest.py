from pathlib import Path
import hashlib, sys
ROOT=Path(__file__).resolve().parents[1]
M=ROOT/'MANIFEST.sha256'
errs=[]; n=0
for line in M.read_text().splitlines():
    if not line.strip(): continue
    want, rel=line.split('  ',1); p=ROOT/rel
    if not p.is_file(): errs.append(f'missing:{rel}'); continue
    h=hashlib.sha256(p.read_bytes()).hexdigest(); n+=1
    if h!=want: errs.append(f'hash:{rel}')
print({'manifest_files_checked':n,'errors':errs,'pass':not errs})
sys.exit(0 if not errs else 1)
