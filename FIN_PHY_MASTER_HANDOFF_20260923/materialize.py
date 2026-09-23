#!/usr/bin/env python3
from pathlib import Path
import argparse, json, shutil, hashlib, sys
ROOT=Path(__file__).resolve().parent

def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def materialize(cid,dest):
 fm=ROOT/'campaigns'/cid/'FILEMAP.json'
 if not fm.exists(): raise SystemExit(f'unknown campaign: {cid}')
 data=json.loads(fm.read_text())
 dest=Path(dest); dest.mkdir(parents=True,exist_ok=True)
 for x in data['files']:
  src=ROOT/x['store_path']; out=dest/x['path']; out.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(src,out)
  if sha(out)!=x['sha256']: raise SystemExit(f'hash mismatch after copy: {cid}/{x["path"]}')
 print(f'{cid}: materialized {len(data["files"])} files -> {dest}')

ap=argparse.ArgumentParser(); g=ap.add_mutually_exclusive_group(required=True)
g.add_argument('--campaign'); g.add_argument('--all',action='store_true'); ap.add_argument('--destination',required=True)
a=ap.parse_args()
if a.all:
 for d in sorted((ROOT/'campaigns').iterdir()):
  if d.is_dir(): materialize(d.name,Path(a.destination)/d.name)
else: materialize(a.campaign,a.destination)
