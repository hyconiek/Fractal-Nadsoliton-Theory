#!/usr/bin/env python3
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parent
checks=(root/'MANIFEST.sha256').read_text().splitlines()
for row in checks:
    sha, filename=row.split('  ',1)
    assert (root/filename).is_file(), filename
    assert hashlib.sha256((root/filename).read_bytes()).hexdigest()==sha,filename
print(f'AH package hashes: {len(checks)}/{len(checks)} PASS')
