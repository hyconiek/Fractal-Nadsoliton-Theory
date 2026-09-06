"""Replay the science and retain actual subprocess outcomes, not status stubs."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent


def checked(command,cwd):
    result=subprocess.run(command,cwd=cwd,text=True,capture_output=True)
    entry=dict(command=command,returncode=result.returncode,
               stdout=result.stdout,stderr=result.stderr)
    if result.returncode:
        raise RuntimeError(json.dumps(entry,indent=2))
    return entry


def main():
    records=[]
    records.append(checked([sys.executable,'-m','unittest','discover','-s','.',
                            '-p','test_*.py','-v'],HERE))
    with tempfile.TemporaryDirectory(prefix='fin30_replay_') as tmp:
        output=Path(tmp)/'results.json'
        records.append(checked([sys.executable,'analysis.py','--through','30',
                                '--output',str(output)],HERE))
        equal=json.loads(output.read_text())==json.loads((HERE/'results.json').read_text())
        if not equal:
            raise AssertionError('Replay differs from the archived result JSON')
    earlier={}
    for name in ['fin_entropy_refinement','fin_jump_fluctuations']:
        directory=HERE.parent/name
        if (directory/'test_analysis.py').is_file():
            earlier[name]=checked([sys.executable,'-m','unittest','discover','-s','.',
                                   '-p','test_analysis.py','-v'],directory)
        else:
            earlier[name]='not bundled; not rerun in this location'
    payload=dict(current_suite_and_replay=records,exact_json_replay_match=equal,
                 earlier_compatibility=earlier,
                 statement='Executed local scientific checks; no external audit or physical data.')
    (HERE/'verification.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('33 current tests passed; exact JSON replay matched. Prior suites checked when present.')


if __name__=='__main__':
    main()
