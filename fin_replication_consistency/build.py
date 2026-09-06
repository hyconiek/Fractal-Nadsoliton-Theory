"""Verify, render and package the completed finite study; no network access."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import zipfile

os.environ.setdefault('MPLCONFIGDIR','/tmp/fin-replication-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

import analysis as a
import certify as c

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def execute(command,cwd=HERE):
    result=subprocess.run(command,cwd=cwd,capture_output=True,text=True)
    record=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    if result.returncode:raise RuntimeError(json.dumps(record,indent=2))
    return record


def verify():
    records=[execute([sys.executable,'-m','unittest','discover','-s','.',
                       '-p','test_*.py','-v'])]
    old=json.loads((HERE/'results.json').read_text())
    exact=json.loads((HERE/'sign_certificate.json').read_text())
    if a.run()!=old:raise AssertionError('Numerical replay differs from archived results')
    if c.run()!=exact:raise AssertionError('Exact sign certificate does not replay')
    prior=ROOT/'fin_adaptive_30'
    if (prior/'test_analysis.py').exists():
        records.append(execute([sys.executable,'-m','unittest','discover','-s','.',
                                 '-p','test_*.py','-v'],prior))
    payload=dict(executed_checks=records,numerical_replay_equal=True,exact_certificate_equal=True,
                 current_tests=25,prior_tests=33 if len(records)>1 else 'not bundled',
                 no_external_evidence=True)
    (HERE/'verification.json').write_text(json.dumps(payload,indent=2)+'\n')


def figure(result,certificate):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
    s=result['strict_row_sum'];M=np.arange(2,1001);delta=1e-4
    one=np.array([a.finite_bound(s,delta,int(n-1)) for n in M])
    two=delta+(s-delta)/(M//2)
    ax=axes[0,0]
    ax.loglog(M,one,label='One-origin bound')
    ax.loglog(M,two,label='Two-origin Gram bound')
    ax.axhline(delta,color='black',ls=':',label='Unlimited-replication bound')
    ax.set(xlabel='Maximum consistent population size M',ylabel='Upper bound on mixed-origin rate',
           title='A. Extra replication adds constraints')
    ax.legend(fontsize=8)
    ax=axes[0,1]
    image=ax.imshow(result['common_pair_matrix'],cmap='viridis',origin='lower')
    fig.colorbar(image,ax=ax,fraction=.046,label='Joint-change rate')
    ax.set(xlabel='Origin j',ylabel='Origin i',title='B. Finite strict activity: zero diagonal, not PSD')
    ax=axes[1,0]
    cases=result['maximal_activity_counterexamples']
    ax.bar([str(x['eta']) for x in cases],[x['transition_00_to_12'] for x in cases],color=['#2874a6','#d68910'])
    ax.set(xlabel='Mixture parameter',ylabel='Rate from (0,0) to (1,2)',
           title='C. Same maximal activity, different events')
    ax=axes[1,1]
    names=['strict_cycle','legacy_cycle','legacy_line'];x=np.arange(3)
    ax.bar(x-.17,[certificate[n]['fixed_points'] for n in names],.34,label='Fixed points')
    ax.bar(x+.17,[certificate[n]['two_cycles'] for n in names],.34,label='Two-cycles')
    ax.set_xticks(x,['Strict cycle','Legacy cycle','Legacy line'])
    ax.set(ylabel='Exact certified count',title='D. Supplied synchronous threshold rule')
    ax.legend()
    fig.savefig(HERE/'findings.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(HERE/'findings.png',dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true')
    p.add_argument('--pdf',action='store_true');p.add_argument('--package',action='store_true')
    opts=p.parse_args()
    result=json.loads((HERE/'results.json').read_text())
    certificate=json.loads((HERE/'sign_certificate.json').read_text())
    if opts.verify:verify()
    figure(result,certificate)
    if opts.pdf:
        logs=[execute(['pdflatex','-interaction=nonstopmode','-halt-on-error','report.tex']) for _ in range(2)]
        (HERE/'pdf_build.json').write_text(json.dumps(logs,indent=2)+'\n')
    inputs=['AGENTS.md','DIAGRAMS_KERNEL_TRANSFORMATION.md',
        'FIN_Programs_31_40_Negative_Information_Coupling.md',
        'FIN_Post_41_50_Methodology_Correction_and_Mirror_Coupling.md',
        'The FIN Kernel as an Unknown Mat.md','neural_coupling_investigation.py',
        'fin_adaptive_30/report.tex']
    hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() if (ROOT/n).is_file()
            else 'not available in standalone replay' for n in inputs}
    provenance=dict(campaign='ST8579--ST8590',date='2026-09-06',
        scope='Selected source-level audit; no full repository revalidation or external evidence.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        input_sha256=hashes,proof_decisions='Exact rational/integer certificate plus analytic proofs; remaining decimals are floating checks.')
    dest='provenance.json' if (ROOT/'AGENTS.md').is_file() else 'replay_provenance.json'
    (HERE/dest).write_text(json.dumps(provenance,indent=2)+'\n')
    paths=[f for f in sorted(HERE.iterdir()) if f.is_file() and f.suffix in
           {'.py','.md','.json','.tex','.pdf','.png'}]
    lines=[hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.name for f in paths]
    (HERE/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')
    if opts.package:
        target=ROOT/'FIN_ST8579_ST8590_Replication_Consistency_Package.zip'
        with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for f in paths+[HERE/'SHA256SUMS.txt']:
                z.write(f,'fin_replication_consistency/'+f.name)
        digest=hashlib.sha256(target.read_bytes()).hexdigest()
        target.with_suffix('.sha256').write_text(digest+'  '+target.name+'\n')
        with zipfile.ZipFile(target) as z:
            if z.testzip() is not None:raise RuntimeError('ZIP integrity failure')
            with tempfile.TemporaryDirectory(prefix='fin_replication_replay_') as tmp:
                z.extractall(tmp)
                work=Path(tmp)/'fin_replication_consistency'
                execute(['sha256sum','-c','SHA256SUMS.txt'],work)
                execute([sys.executable,'-m','unittest','discover','-s','.',
                         '-p','test_*.py'],work)
                check=execute([sys.executable,'analysis.py','--output','replay.json'],work)
                if json.loads((work/'replay.json').read_text())!=result:
                    raise AssertionError('Clean package numerical replay mismatch')
        print(target)
        print('Verified clean extraction, 25 scientific tests and exact numerical replay.')
    print('Built report support and manifest.')


if __name__=='__main__':main()
