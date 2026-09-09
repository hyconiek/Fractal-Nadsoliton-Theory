"""Replay the finite discovery package without generating a PDF."""
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
    from . import research,operational
except ImportError:
    import research,operational


HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def source_tex_check(path):
    source=path.read_text()
    cleaned=re.sub(r'(?<!\\)%[^\n]*','',source)
    literal_removed=re.sub(r'\\[{}]','',cleaned)
    depth=0
    for char in literal_removed:
        if char=='{':depth+=1
        elif char=='}':depth-=1
        if depth<0:raise AssertionError('Unbalanced closing TeX brace')
    assert depth==0
    stack=[]
    for action,env in re.findall(r'\\(begin|end)\{([^}]+)\}',cleaned):
        if action=='begin':stack.append(env)
        else:assert stack and stack.pop()==env
    assert not stack
    labels=set(re.findall(r'\\label\{([^}]+)\}',source))
    refs=set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',source))
    assert refs<=labels
    bib=set(re.findall(r'\\bibitem\{([^}]+)\}',source))
    cites={item for group in re.findall(r'\\cite\{([^}]+)\}',source) for item in group.split(',')}
    assert cites<=bib
    return dict(braces_balanced=True,environments_balanced=True,references_resolved=True,
        citations_resolved=True,compilation_performed=False)


def main():
    run=subprocess.run([sys.executable,'-m','unittest',
        'fin_chiral_selection.test_research','fin_chiral_selection.test_operational','-v'],
        cwd=ROOT,text=True,capture_output=True)
    output=run.stdout+run.stderr
    if run.returncode:raise RuntimeError(output)
    match=re.search(r'Ran (\d+) tests',output)
    assert match and int(match.group(1))==17
    comparisons={}
    for filename,producer in [('results.json',research.run),('operational_results.json',operational.run)]:
        saved=json.loads((HERE/filename).read_text())
        current=producer()
        assert current==saved,filename+' does not replay exactly'
        comparisons[filename]=True
    report=ROOT/'FIN_Adaptive_Feedback_Operational_Obstruction.tex'
    tex=source_tex_check(report)
    assert not report.with_suffix('.pdf').exists()
    assert not list(HERE.rglob('*.pdf'))
    sources=[report,ROOT/'The FIN Kernel as an Unknown Mat.md',
             ROOT/'fin_replication_consistency/certify.py',
             ROOT/'fin_projected_learning/research.py',HERE/'PROOF.md',HERE/'WORKLOG.md']
    sources+=list(HERE.glob('*.py'))
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    audit=dict(status='PASS',scientific_tests=17,test_output=output,
        exact_result_replay=comparisons,tex_static_checks=tex,no_pdf_generated=True,
        versions=dict(python=platform.python_version(),numpy=numpy.__version__,
            scipy=scipy.__version__,sympy=sympy.__version__),source_sha256=hashes,
        completion_audit=dict(
            actual_update_source='Original pure equation inspected; mixed extension kept separate.',
            new_important_result='Exact finite operational non-affinity obstruction, including strict positive edges and canonical vertex readout.',
            falsification='Zero learning, insensitive ensembles, Laplacian variant, positive cone, exact two-site independent solution and primary literature critique checked.',
            prior_repository_scope='Builds on prior broad state/source/composition audits; no claim to reprove every archived result.',
            legacy_strict_split='Arbitrary-K0 coefficient applies separately; no kernel substitution or role transfer.',
            physical_scope='Standard mixture/branchwise readout premises explicit; mean-field and nonstandard alternatives remain open.',
            final_report=str(report.relative_to(ROOT)),
            final_frontier='Sourced finite microscopic composition/controller law and a controlled operational/mean-field limit.'),
        limitations=['No laboratory evidence or dimensional calibration.',
            'No source derivation, selector, legacy bridge, role transfer, SM/GR or ToE closure.',
            'No global priority claimed for Gisin-type nonlinearity or saddle-focus winding.',
            'No proof of the fixed-K0 fixed-phase logarithmic asymptotic or full-system chiral stability.'])
    (HERE/'verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_sha256','test_output']},indent=2))


if __name__=='__main__':main()
