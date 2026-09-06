"""Build figures, the adaptive ledger and a reproducibility manifest.

Does not certify theorems. Run after analysis.py and before LaTeX. Run once
more after the PDF is built to update hashes. No network or package install.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

os.environ.setdefault('MPLCONFIGDIR','/tmp/fin-adaptive30-mpl')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.linalg import eigh,expm

import analysis as a

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def latex(s):
    replacements={"\\":r"\textbackslash{}", "&":r"\&", "%":r"\%",
                  "$":r"\$", "#":r"\#", "_":r"\_", "{":r"\{",
                  "}":r"\}", "~":r"\textasciitilde{}", "^":r"\textasciicircum{}",
                  "<":r"\textless{}", ">":r"\textgreater{}"}
    return ''.join(replacements.get(c,c) for c in s)


def figures():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    W,V=a.kernels(); A=a.lap(W)
    fig,axes=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
    d=np.arange(1,7)
    ax=axes[0,0]
    ax.plot(d,W[0,1:7],'o-',label='Strict W')
    ax.plot(d,V[0,1:7],'s-',label='Canonical legacy V')
    ax.axhline(0,color='black',lw=.7)
    ax.set(xlabel='Cyclic distance',ylabel='Signed weight',title='A. Different kernels; no silent replacement')
    ax.legend()
    ts=np.geomspace(.0001,.2,80)
    ax=axes[0,1]
    heat=np.array([expm(-t*A)[0,1] for t in ts])
    born=np.array([abs(expm(-1j*t*A)[0,1])**2 for t in ts])
    ax.loglog(ts,heat,label='Heat vertex probability')
    ax.loglog(ts,born,label='Unitary Born vertex probability')
    ax.set(xlabel='Dimensionless t',ylabel='Transition probability',title='B. One generator, different record laws')
    ax.legend()
    ax=axes[1,0]
    ts=np.geomspace(.0001,.3,100)
    for theta in [0.,.35]:
        G,xs=a.population_generator(W,theta=theta)
        vals,vecs=eigh(G)
        i,j=xs.index((0,0)),xs.index((1,1))
        probs=np.exp(np.outer(ts,vals))@(vecs[i]*vecs[j])
        ax.loglog(ts,probs,label=f'Common-event fraction {theta:g}')
    ax.set(xlabel='Dimensionless t',ylabel='P[(0,0) to (1,1)]',title='C. Identical singleton laws, distinct pair laws')
    ax.legend()
    ax=axes[1,1]
    deltas=np.geomspace(1e-4,.08,100)
    q=np.array([a.observed_double_rate(.7,.2,.6,t)[2] for t in deltas])
    ax.loglog(deltas,.5*(1.4*deltas)**2,label='Conservative single-change null bound')
    ax.loglog(deltas,q,label='Resolved common-flip alternative')
    ax.axvline(.01,color='gray',ls=':',lw=1)
    ax.set(xlabel='Bin width in model clock',ylabel='P[at least two detected changes]',
           title='D. Finite-resolution test (supplied detector)')
    ax.legend(fontsize=8)
    fig.savefig(HERE/'findings.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(HERE/'findings.png',dpi=160)
    plt.close(fig)


def main():
    rows=json.loads((HERE/'results.json').read_text())
    if [r['round'] for r in rows] != list(range(1,31)):
        raise ValueError('All thirty real result records are required for the final report')
    tex=[]; md=['# FIN: Thirty-step research ledger','']
    for r in rows:
        tex.append(r'\subsection*{'+latex(f"{r['round']:02d}. {r['program']} -- {r['question']}")+'}')
        tex.append(r'\noindent\textbf{Outcome.} '+latex(r['result'])+r'\par')
        tex.append(r'\noindent\textbf{Status.} '+latex(r['status'])+r'\par')
        tex.append(r'\noindent\textbf{Next-step rationale.} '+latex(r['next_reason'])+r'\par\medskip')
        md += [f"## {r['round']:02d}. {r['program']} — {r['question']}",'',
               r['result'],'',f"Status: {r['status']}",'',
               f"Next-step rationale: {r['next_reason']}",'',
               'Numerical/check evidence:','', '```json',json.dumps(r['evidence'],indent=2),'```','']
    (HERE/'rounds.tex').write_text('\n'.join(tex)+'\n')
    (HERE/'RESEARCH_LEDGER_EN.md').write_text('\n'.join(md))
    figures()
    source_names=['AGENTS.md','SUMMARY_GROK.md',
        'FIN_Theory_Compendium_From_Fractal_Information_to_Current_Mathematics_EN.tex',
        'fin_entropy_refinement/PROOFS_EN.md','fin_entropy_refinement/analysis.py',
        'fin_jump_fluctuations/PROOFS_EN.md','fin_jump_fluctuations/analysis.py']
    sources={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() if (ROOT/s).is_file()
             else 'unavailable in this replay location' for s in source_names}
    try:
        files=subprocess.check_output(['rg','--files','-g','!hermes-agent/**','-g','!node_modules/**',
                                       '-g','!.git/**'],cwd=ROOT,text=True).splitlines()
    except (OSError,subprocess.CalledProcessError):
        files=[p.name for p in HERE.iterdir() if p.is_file()]
    ext={}
    for f in files:
        suffix=Path(f).suffix or '<none>'
        ext[suffix]=ext.get(suffix,0)+1
    try:
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True,
                                     stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError):
        head='unavailable in this replay location'
    provenance=dict(campaign='ST8549--ST8578',rounds=30,random_seed=8549,
        git_head=head,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                         matplotlib=matplotlib.__version__),
        initial_inventory_count_approx=33000,current_inventory_by_extension=ext,
        scope='Repository-wide file/section inventory; selected proof-level reading, not revalidation of all historical claims.',
        input_sha256=sources,
        local_commands=['python3 analysis.py --through 30',
            "python3 -m unittest discover -s . -p 'test_*.py' -v",'python3 build_artifacts.py',
            'pdflatex -interaction=nonstopmode -halt-on-error report.tex (twice)'],
        precision='Analytic proofs with rational/symbolic checks; NumPy/SciPy floating computations are not interval certificates.',
        no_physical_evidence=True)
    # Preserve the supplied source manifest when replaying a standalone bundle.
    destination='provenance.json' if (ROOT/'AGENTS.md').is_file() else 'replay_provenance.json'
    (HERE/destination).write_text(json.dumps(provenance,indent=2)+'\n')
    wanted={'.py','.json','.md','.tex','.pdf','.png'}
    hashes=[]
    for path in sorted(HERE.iterdir()):
        if path.is_file() and path.suffix in wanted:
            hashes.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+path.name)
    (HERE/'SHA256SUMS.txt').write_text('\n'.join(hashes)+'\n')
    print(f'Built thirty-round ledger, figure and manifest ({len(hashes)} artifact hashes).')


if __name__=='__main__':
    main()
