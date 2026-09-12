"""Final evidence audit of the new quantum-lift research goal, without PDF."""
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
    from . import research,collisions
except ImportError:
    import research,collisions


HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def check_tex(path):
    text=path.read_text()
    cleaned=re.sub(r'(?<!\\)%[^\n]*','',text)
    cleaned=re.sub(r'\\[{}]','',cleaned)
    depth=0
    for char in cleaned:
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
    references=set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text))
    assert references<=set(labels)
    return dict(braces=True,environments=True,references=True,compilation_performed=False)


def main():
    result=subprocess.run([sys.executable,'-m','unittest','fin_quantum_lift.test_research',
        'fin_quantum_lift.test_collisions','-v'],cwd=ROOT,text=True,capture_output=True)
    log=result.stdout+result.stderr
    if result.returncode:raise RuntimeError(log)
    match=re.search(r'Ran (\d+) tests',log)
    assert match and int(match.group(1))==25
    replay={}
    for filename,producer in [('results.json',research.run),('collision_results.json',collisions.run)]:
        assert json.loads((HERE/filename).read_text())==producer(),filename
        replay[filename]='exact replay'
    report=ROOT/'FIN_Quantum_Source_Correlation_and_Programming_Report.tex'
    tex=check_tex(report)
    assert not report.with_suffix('.pdf').exists()
    assert not list(HERE.rglob('*.pdf'))
    sources=list(HERE.glob('*.py'))+[HERE/'PROOF.md',HERE/'WORKLOG.md',report,
        ROOT/'The FIN Kernel as an Unknown Mat.md',
        ROOT/'fin_projected_learning/research.py',ROOT/'fin_projected_learning/geometry.py',
        ROOT/'fin_replication_consistency/certify.py',
        ROOT/'FIN_Programs_507_516_Nonlinear_Source_and_Certified_Localization_Report_EN.tex']
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(sources)}
    output=dict(status='PASS',scientific_tests=25,test_output=log,result_replay=replay,
        tex_static_checks=tex,no_pdf_generated=True,
        versions=dict(python=platform.python_version(),numpy=numpy.__version__,
            scipy=scipy.__version__,sympy=sympy.__version__),source_sha256=hashes,
        completion_audit=dict(
            original_objective='Autonomous rigorous FIN research until an important finding, then one report; not another fixed-round campaign.',
            repository_and_guardrails='Original projected source, preceding operational/state/composition audits, exact strict enclosures and distinct canonical legacy reference consulted.',
            new_result='Trivial two-local commutant implies maximally mixed full-rank product stationarity on every finite connected graph with canonical V edges, even with arbitrary local fields.',
            significance='It rules out a whole independent microscopic equilibrium class, not merely one numerical seed; stationary correlated storage is separately shown insufficient for propagation.',
            falsification='n=2 and rank-one product exceptions; local-field cancellation check; invalid correlated completion including optimal-loading boundary; flat-band pulse; hidden program noise; variance-only false positive; signed legacy and arithmetic-scope split.',
            positive_structure='Exact two-copy quantum dynamics, explicit correlated marginals, short-time BBGKY bound and a reference-assisted finite-copy CPTP simulation bound.',
            literature='Primary mean-field, sample-based simulation and no-programming sources inspected; their known results are not counted as new discoveries.',
            physical_scope='Tensor law, state/correlation supply, fresh copies and clock are stated inputs. No claim that a compatible simulation derives a fundamental FIN law.',
            final_report=str(report.relative_to(ROOT)),
            next_direction='A sourced microscopic law and correlation/program preparation that predicts both propagation and finite operational corrections.'),
        limitations=['General analytical proofs are not proof-assistant formalizations.',
            'The full-rank theorem does not exclude rank-deficient product or correlated stationary states.',
            'Programmed approximation is not exact finite-eta self-feedback or source-independent kernel emergence.',
            'Bounds and loading optimization concern the declared processor, not all quantum algorithms.',
            'No laboratory, units, selector, legacy bridge, physical-role, SM/GR or ToE closure.'])
    (HERE/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['source_sha256','test_output']},indent=2))


if __name__=='__main__':main()
