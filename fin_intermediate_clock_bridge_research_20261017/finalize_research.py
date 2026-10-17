#!/usr/bin/env python3
"""New packaging helper only; preserve all comparison scientific sources."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    result=json.loads((HERE/'RESULTS.json').read_text())
    assert all(c['passed'] for c in result['checks']) and len(result['checks'])==12
    for path in HERE.glob('*.py'):ast.parse(path.read_text())
    for path in HERE.glob('*.json'):json.loads(path.read_text(),parse_constant=lambda v:(_ for _ in ()).throw(ValueError(v)))
    for path in HERE.glob('*.md'):
        for t in re.findall(r'\[[^\]]+\]\(([^)]+)\)',path.read_text()):
            if not t.startswith(('http:','https:','#')):assert (path.parent/t).exists(),(path,t)
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    print(json.dumps({'status':'PASS_COMPARISON_APPENDIX_PACKAGE','entries':len(files),'read_only':a.check}))
if __name__=='__main__':main()
