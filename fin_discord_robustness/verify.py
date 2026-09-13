"""Final source-only replay and scope audit for the active discovery study."""
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys

import numpy
import scipy
import sympy

try:
    from . import research
except ImportError:
    import research

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def check_tex(path):
    text=path.read_text();clean=re.sub(r'(?<!\\)%[^\n]*','',text)
    clean=re.sub(r'\\[{}]','',clean);depth=0
    for char in clean:
        if char=='{':depth+=1
        elif char=='}':depth-=1
        assert depth>=0
    assert depth==0
    stack=[]
    for action,name in re.findall(r'\\(begin|end)\{([^}]+)\}',text):
        if action=='begin':stack.append(name)
        else:assert stack and stack.pop()==name
    assert not stack
    labels=re.findall(r'\\label\{([^}]+)\}',text)
    assert len(labels)==len(set(labels))
    assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text))<=set(labels)
    return dict(braces=True,environments=True,references=True,compilation_performed=False)


def main():
    result=subprocess.run([sys.executable,'-m','unittest','fin_discord_robustness.test_research','-v'],
        cwd=ROOT,text=True,capture_output=True)
    log=result.stdout+result.stderr
    if result.returncode:raise RuntimeError(log)
    assert re.search(r'Ran 20 tests',log)
    assert research.run()==json.loads((HERE/'results.json').read_text())
    report=ROOT/'FIN_Discord_Robustness_and_Operational_Identifiability.tex'
    tex=check_tex(report)
    assert not report.with_suffix('.pdf').exists() and not list(HERE.rglob('*.pdf'))
    sources=list(HERE.glob('*.py'))+[HERE/'PROOF.md',HERE/'WORKLOG.md',report,
        ROOT/'fin_quantum_lift/research.py',ROOT/'fin_projected_learning/research.py',
        ROOT/'fin_replication_consistency/certify.py']
    data=dict(status='PASS',scientific_tests=20,test_output=log,results_exactly_replayed=True,
        tex_static_checks=tex,no_pdf_generated=True,
        versions=dict(python=platform.python_version(),numpy=numpy.__version__,
                      scipy=scipy.__version__,sympy=sympy.__version__),
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
        completion_audit=dict(
            objective='New rigorous research goal after the completed separability/necessary-discord report.',
            frontier='The explicit unresolved microscopic exchange shift, then the operational and preparation ambiguities it exposed.',
            source_and_guards='Full-source, mixed-flow and pure-flow premises remain separate; strict and legacy programs are checked without role transfer.',
            important_result='Universal stationary preparation channels have identical complete local channels but are separable exactly at balanced mixing and NPT otherwise; local FIN data cannot identify this joint resource.',
            additional_proofs='All mixed-flow-equivalent Hamiltonians retain the non-CQ bound; the pure-flow family obeys an asymmetry/discord/stationarity tradeoff.',
            falsification='True spectral transfer exceptions, CQ destroyed by state averaging, unequal versus equal marginals, controlled-U versus passive access, non-LOCC universal routing, zero-kernel entanglement, and independent legacy program.',
            evidence='Exact polynomial identities, symbolic CPTP marginal formulas, exact negative partial-transpose minors and rational bounds; numerical witnesses separately replayed.',
            literature='Primary controlled-operation and cloning sources used to bound interpretation; no global novelty or universal optimal cloning claim.',
            final_report=str(report.relative_to(ROOT)),
            physical_nonclosure='Program, joint preparation law, factor labels, control/energy access and calibration are not derived from FIN.',
            next_direction='A sourced joint preparation channel and observable/control algebra, with independent joint and coherent-control predictions.'),
        limitations=['The asymmetric CQ witness matches only the average marginal, not two individual strict marginals.',
            'Passive operational equivalence excludes controlled-Hamiltonian and energy-query access.',
            'The separable one-copy channel is global; universal fixed-party LOCC requires additional quantum resources.',
            'Stationary storage does not derive the strict propagator or a physical source.',
            'No selector, legacy bridge, dimensional unit, laboratory, SM/GR or ToE closure.'])
    (HERE/'verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['source_sha256','test_output']},indent=2))


if __name__=='__main__':main()
