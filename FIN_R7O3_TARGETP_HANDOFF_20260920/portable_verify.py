from pathlib import Path
import subprocess, sys
ROOT=Path(__file__).resolve().parent
cmds=[
 [sys.executable,'verify.py'],
 [sys.executable,'src/global_join_checker.py'],
 [sys.executable,'src/hostile_mutations.py'],
 [sys.executable,'src/verify_manifest.py'],
]
for cmd in cmds:
    print('RUN', ' '.join(cmd), flush=True)
    p=subprocess.run(cmd,cwd=ROOT)
    if p.returncode:
        print('portable_replay=FAIL',cmd,p.returncode); raise SystemExit(p.returncode)
print('portable_replay=PASS')
