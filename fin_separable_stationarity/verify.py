"""Replay the separable/discord proof package without generating a PDF."""
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
    from . import research,discord
except ImportError:
    import research,discord

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def tex_check(path):
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
    result=subprocess.run([sys.executable,'-m','unittest','fin_separable_stationarity.test_research','-v'],
        cwd=ROOT,text=True,capture_output=True)
    log=result.stdout+result.stderr
    if result.returncode:raise RuntimeError(log)
    assert re.search(r'Ran 18 tests',log)
    replay={}
    for name,producer in [('results.json',research.run),('discord_results.json',discord.run)]:
        assert producer()==json.loads((HERE/name).read_text()),name
        replay[name]='exact replay'
    report=ROOT/'FIN_Separable_Stationarity_and_Necessary_Discord.tex'
    checks=tex_check(report)
    assert not report.with_suffix('.pdf').exists() and not list(HERE.rglob('*.pdf'))
    sources=list(HERE.glob('*.py'))+[HERE/'PROOF.md',HERE/'WORKLOG.md',report,
        ROOT/'fin_quantum_lift/research.py',ROOT/'fin_projected_learning/research.py',
        ROOT/'fin_replication_consistency/certify.py']
    data=dict(status='PASS',scientific_tests=18,test_output=log,result_replay=replay,
        tex_static_checks=checks,no_pdf_generated=True,
        versions=dict(python=platform.python_version(),numpy=numpy.__version__,
                      scipy=scipy.__version__,sympy=sympy.__version__),
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
        completion_audit=dict(
            objective='New rigorous FIN research goal after the completed correlation-floor study.',
            source_and_guards='Canonical interaction, exact strict inputs and a separately checked canonical legacy reference; no source/role transfer.',
            important_result='Explicit separable stationary strict state exists, but every canonical stationary state with the declared real strict marginal is non-CQ on that side.',
            positive_mechanism='22 product-state decomposition and explicit local projective success instrument with probability 1/2, using an amplified supplied program.',
            mathematical_evidence='Integer Hadamard/cut identities, rational positivity/gap bounds, exact block inverse, general no-CQ proof and trace-distance bounds.',
            falsification='Zero loading, non-Hermitian/overloaded program rejection, preserved record, PPT versus separability, zero commutator versus discord, SWAP alternative, exceptional exchange shift, and inequivalent fourth-moment frames.',
            literature='Known discord and lazy-state results are cited as tools and interpretation boundaries, not counted as new concepts.',
            nonclosure='Program, cut law, heralding and physical identification remain inputs; no autonomous FIN source or laboratory evidence is derived.',
            report=str(report.relative_to(ROOT)),
            next_direction='An independently justified amplified coherent program and higher-order preparation law with testable joint predictions.'),
        limitations=['No entropic discord value or universal entanglement requirement is claimed.',
            'No-CQ theorem requires the stated invertible block and simple unequal marginal eigenvalues.',
            'Strict quantitative bounds do not automatically cover all Hartree-equivalent interactions or legacy.',
            'The cut threshold is sufficient, and cut-count optimality is architecture scoped.',
            'No physical clock, apparatus, selector, bridge, role transfer, SM/GR or ToE closure.'])
    (HERE/'verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['test_output','source_sha256']},indent=2))


if __name__=='__main__':main()
