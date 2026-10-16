#!/usr/bin/env python3
"""MM package file integrity checker only, not scientific proof."""
import hashlib
from pathlib import Path
R=Path(__file__).resolve().parent
exclude={'MANIFEST.sha256','PACKAGE_VERIFICATION.json','verify_package.py','build_mm_documents.py'}
lines=(R/'MANIFEST.sha256').read_text().splitlines()
expected={}
for line in lines:
    h,relative=line.split('  ',1)
    if relative in expected:raise SystemExit(f'Duplicate path: {relative}')
    expected[relative]=h
real={str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and p.name not in exclude and '__pycache__' not in p.parts}
if real != set(expected):
    raise SystemExit(f'Coverage ERROR missing={sorted(real-set(expected))}; surplus={sorted(set(expected)-real)}')
for relative,h in expected.items():
    actual=hashlib.sha256((R/relative).read_bytes()).hexdigest()
    if actual != h: raise SystemExit(f'HASH ERROR: {relative}')
print(f'MANIFEST OK {len(expected)}/{len(expected)}; coverage complete for declared scientific source files.')
