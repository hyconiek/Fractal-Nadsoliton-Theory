#!/usr/bin/env python3
from __future__ import annotations
import subprocess
from pathlib import Path

FAR=Path(__file__).resolve().parent
GEN=FAR/"generated"

def main()->None:
    targets=sorted([p for p in GEN.glob("p13*.json")]+[p for p in GEN.glob("p13*.csv")])
    if not targets:
        raise SystemExit("No p13* artifacts found")
    missing=[]
    for t in targets:
        rel=t.as_posix()
        r=subprocess.run(["git","ls-files","--error-unmatch",rel],capture_output=True,text=True)
        if r.returncode!=0:
            missing.append(rel)
    if missing:
        raise SystemExit("Untracked artifact(s):\n"+"\n".join(missing))
    print(f"[P1365] tracked artifacts: {len(targets)}")

if __name__=="__main__":
    main()
