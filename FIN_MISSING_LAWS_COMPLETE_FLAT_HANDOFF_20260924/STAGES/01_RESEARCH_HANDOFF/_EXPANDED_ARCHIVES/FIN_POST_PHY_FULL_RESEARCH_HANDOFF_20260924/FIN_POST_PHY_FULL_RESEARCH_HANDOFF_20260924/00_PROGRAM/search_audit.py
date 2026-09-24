"""Repo-wide content de-dup index for the post-PHY design; no research producers."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1]; HERE=Path(__file__).resolve().parent
MASTER=ROOT/'FIN_PHY_MASTER_HANDOFF_20260923'
GROUPS={
 'proof_source':'theorem|no[- ]go|source|selector|rank|carrier|failed approaches|next (tasks|atoms|steps)',
 'refinement_geometry':'refinement|continuum|locality|geometry|Dirichlet|spectral|resonance|symmetry|phase',
 'memory_inference':'Hankel|realization|hidden.memory|HMM|minimax|e-process|anytime|identifiability|calibration',
 'dynamics_physics':'causal|propagator|dynamics|Markov|unitary|wave|quantum|gravity|measurement|falsification|operational|physical.bridge',
 'new_bridge_candidates':'Dirichlet.to.Neumann|two.port|impedance|nonnegative.rank|positive.realization|predictable.*calibr|adaptive.*Hankel|convex.hull.*null|spectral.factorization|near.resonance|phase.slip|phase.locking.*(gap|uniform)|Mosco|random.geometric|RNG|automorphism.*perturb|memory.*refinement'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out={'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
         'agents_sha256':sha(ROOT/'AGENTS.md'),'method':'rg content search over Markdown/Python/JSON plus separately indexed content-store logical files; default ignored/binary files excluded','groups':{}}
    for name,pattern in GROUPS.items():
        args=['rg','-l','-i',pattern,'-g','*.md','-g','*.py','-g','*.json','-g','!fin_post_phy_review/**','-g','!FIN_POST_PHY*','.']
        p=subprocess.run(args,cwd=ROOT,text=True,capture_output=True,check=False)
        assert p.returncode in (0,1),p.stderr
        files=sorted(p.stdout.splitlines())
        out['groups'][name]={'pattern':pattern,'file_count':len(files),'files':files}
        print(name,len(files),flush=True)
    mapped={}
    for fm in sorted((MASTER/'campaigns').glob('*/FILEMAP.json')):
        for r in json.loads(fm.read_text())['files']:
            if Path(r['path']).suffix in ['.md','.py','.json']:
                mapped[fm.parent.name+'::'+r['path']]={'sha256':r['sha256'],'path':str((MASTER/r['store_path']).relative_to(ROOT))}
    out['logical_sources']=mapped
    (HERE/'search_index.json').write_text(json.dumps(out,indent=2)+'\n')
    print('logical text sources',len(mapped))
if __name__=='__main__':main()
