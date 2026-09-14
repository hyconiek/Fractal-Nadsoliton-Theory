"""Rebuild audit outputs and run new plus relevant inherited regressions."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from fin_handoff_audit import research


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--record',action='store_true',help='Write generated audit ledgers after verification.')
    args=parser.parse_args()
    suites=['fin_handoff_audit.test_research','fin_hartree_equivalence.test_research',
            'fin_separable_stationarity.test_research','fin_discord_robustness.test_research']
    logs={}
    for suite,expected in zip(suites,[19,18,18,20]):
        proc=subprocess.run([sys.executable,'-m','unittest',suite,'-v'],cwd=HERE.parent,
                            text=True,capture_output=True,timeout=180)
        logs[suite]=proc.stdout+proc.stderr
        if proc.returncode: raise RuntimeError(logs[suite])
        assert re.search(rf'Ran {expected} tests\b',logs[suite])
        print(suite+': PASS',flush=True)
    results=research.run()
    # Exact layer must replay exactly. Float layer is separately regression-tested.
    if not args.record:
        saved=json.loads((HERE/'results.json').read_text())
        assert results['exact']==saved['exact']
        assert results['manifest']==saved['manifest']
    paths=list(HERE.glob('*.py'))+[HERE/'REPORT.md',HERE/'PROOF.md',
        research.BUNDLE/'FIN_full_chat_research_handoff_pre_and_post_Discord.md',
        research.BUNDLE/'MANIFEST.txt',HERE.parent/'fin_projected_learning/research.py',
        HERE.parent/'fin_replication_consistency/certify.py']
    for package in ['fin_hartree_equivalence','fin_separable_stationarity','fin_discord_robustness']:
        paths.extend((HERE.parent/package).glob('*.py'))
    summary={'status':'PASS','new_tests':19,'inherited_regression_tests':56,
             'python':platform.python_version(),'test_logs':logs,
             'sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in paths},
             'limits':['No replay of absent historical global cover or absent 60-root catalog.',
                       'Exact certificates cover the stated face only.',
                       'The full-seven-coordinate everywhere index bound is refuted, not the stationary-point census.',
                       'No physical/source/selector/legacy-role promotion.']}
    if args.record:
        (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        (HERE/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['sha256','test_logs']},indent=2))


if __name__=='__main__':main()
