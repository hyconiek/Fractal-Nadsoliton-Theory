"""Final replay and requirement audit; no PDF or formal-proof claim."""
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
    tests=subprocess.run([sys.executable,'-m','unittest','fin_hartree_equivalence.test_research','-v'],
                         cwd=ROOT,text=True,capture_output=True)
    log=tests.stdout+tests.stderr
    if tests.returncode:raise RuntimeError(log)
    assert re.search(r'Ran 18 tests',log)
    assert research.run()==json.loads((HERE/'results.json').read_text())
    report=ROOT/'FIN_Hartree_Equivalence_and_Strict_Correlation_Floor.tex'
    tex=check_tex(report)
    assert not report.with_suffix('.pdf').exists() and not list(HERE.rglob('*.pdf'))
    inputs=list(HERE.glob('*.py'))+[HERE/'PROOF.md',HERE/'WORKLOG.md',report,
        ROOT/'fin_quantum_lift/research.py',ROOT/'fin_projected_learning/research.py',
        ROOT/'fin_replication_consistency/certify.py']
    data=dict(status='PASS',scientific_tests=18,test_output=log,results_exactly_replayed=True,
        tex_static_checks=tex,no_pdf_generated=True,
        versions=dict(python=platform.python_version(),numpy=numpy.__version__,
                      scipy=scipy.__version__,sympy=sympy.__version__),
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(inputs)},
        completion_audit=dict(
            objective='New rigorous discovery goal after the completed quantum-source report; previous findings not counted again.',
            source_and_guards='Fixed normal-ordering source, current microscopic frontier, separate strict/legacy references and exact strict spectral providers retained.',
            new_important_result='A rank-free strict commutator witness and a uniform positive correlation floor across all 4357 pure-projector-flow-invisible reciprocal interaction parameters.',
            proof_scope='Exact all-source/mixed/pure classification; analytic rank-free and matching obstructions; exact rational floor backed by an explicit dual witness.',
            falsification='Exceptional exchange shift, bipartite/odd-cycle contrast, actual matching product, mixed-state distinguishability, n=2, zero loading, rank-six chiral states, graph local-field cancellation and non-vacuous stationary correlated state.',
            literature='Commuting-map theory and Pinsker inequality explicitly cited as known tools, not claimed as discoveries.',
            source_nonclosure='Correlation floor is conditional on the stated state/interaction class; no physical preparation, clock, selector, bridge, role transfer or ToE closure.',
            final_report=str(report.relative_to(ROOT)),
            next_direction='A sourced zero-marginal correlation chi satisfying positivity and the noncancellable witness equation, with an operational propagation prediction.'),
        limitations=['Pure projector dynamics is weaker than the full controller source and does not fix controlled phase.',
            'The numerical correlation floor is a two-body result without arbitrary extra local fields.',
            'The floor is not an entanglement witness or an existence theorem for every interaction.',
            'General proofs are not proof-assistant formalizations; finite exact and numerical checks have distinct roles.'])
    (HERE/'verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['source_sha256','test_output']},indent=2))


if __name__=='__main__':main()
