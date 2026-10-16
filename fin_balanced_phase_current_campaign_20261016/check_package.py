#!/usr/bin/env python3
import hashlib, pathlib, sys, zipfile
zpath=pathlib.Path(sys.argv[1]); assert zpath.is_file()
with zipfile.ZipFile(zpath) as z:
 assert z.testzip() is None
 names=set(z.namelist()); root='fin_balanced_phase_current_campaign_20261016/'
 lines=z.read(root+'MANIFEST.sha256').decode().splitlines()
 for row in lines:
  digest,name=row.split('  ',1)
  assert root+name in names, name
  assert hashlib.sha256(z.read(root+name)).hexdigest()==digest,name
 assert len(lines)>=20
 print('ZIP PASS:',len(lines),'scientific output manifest hashes; entries',len(names))
