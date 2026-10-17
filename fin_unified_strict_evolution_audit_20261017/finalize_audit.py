#!/usr/bin/env python3
"""Package only this audit; --check is read-only, no predecessor edits."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
REQUIRED=['STRICT_EVOLUTION_SOURCE_AUDIT.md','AUTONOMOUS_COMPOSITION_THEOREM_OR_OBSTRUCTION.md',
          'PERSISTENT_STRUCTURES_AND_PROPAGATION.md','ENERGY_TRANSFER_AND_CONSERVATION.md',
          'MULTISCALE_FRACTAL_DYNAMICS.md','JOINT_PHENOMENA_EMERGENCE.md',
          'PROOF_GAPS_AND_COUNTEREXAMPLES.md','CLAIM_REGISTER.json','REPORT.md','HANDOFF.md',
          'COMPLETION_AUDIT.md','MODEL_AND_SOURCE_CONTRACT.json','KERNEL_SENSITIVE_TRANSVERSE_SEED.md',
          'LITERATURE_DATABASE.md','LITERATURE_REGISTER.json','BROWSER_READ_LOG.jsonl',
          'REPLAY.md','verify_evolution.py','verification.json','browser_capture.py','finalize_audit.py']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args()
    assert all((HERE/n).is_file() for n in REQUIRED)
    saved=json.loads((HERE/'verification.json').read_text())
    assert saved['status']=='PASS_BOUNDED_SOURCE_AND_EVOLUTION_AUDIT'
    assert len(saved['checks'])==44 and all(v['passed'] for v in saved['checks'])
    assert saved['input_hashes']==saved['source_hashes_after']
    for rel,h in saved['input_hashes'].items():assert sha(ROOT/rel)==h,rel
    assert sha(ROOT/'AGENTS.md')==json.loads((HERE/'MODEL_AND_SOURCE_CONTRACT.json').read_text())['immutable_AGENTS_sha256']
    for p in HERE.glob('*.json'):json.loads(p.read_text(),parse_constant=lambda v:(_ for _ in ()).throw(ValueError(v)))
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    for claim in json.loads((HERE/'CLAIM_REGISTER.json').read_text())['claims']:assert (HERE/claim['evidence']).is_file()
    for p in [*HERE.glob('*.md'),ROOT/'FIN_UNIFIED_EVOLUTION_SOURCE_AND_COMPOSITION_NEXT_CHAIN_20261017.md']:
        assert all(l==l.rstrip() for l in p.read_text().splitlines()),p
        for t in re.findall(r'\[[^\]]+\]\(([^)]+)\)',p.read_text()):
            if not t.startswith(('https:','http:','#')):assert (p.parent/t).exists(),(p,t)
    rows=[json.loads(l) for l in (HERE/'BROWSER_READ_LOG.jsonl').read_text().splitlines()]
    assert {'kink-dynamics-search','kmm-primary-abstract'}<={r['label'] for r in rows}
    for r in rows:
        cache=ROOT/r['optional_cache']
        if cache.exists():assert sha(cache)==r['text_sha256']
    literature=json.loads((HERE/'LITERATURE_REGISTER.json').read_text())
    for r in literature['works']:
        for pkey,hkey in (('temporary_pdf','pdf_sha256'),('temporary_extraction','pdftotext_layout_sha256')):
            p=Path(r[pkey])
            if p.exists():assert sha(p)==r[hkey]
    agents=(ROOT/'AGENTS.md').read_text().splitlines()
    headings=[{'line':i,'heading':l} for i,l in enumerate(agents,1) if l.startswith('## ')]
    index={'date':'2026-10-17','source_hashes':saved['input_hashes'],
           'AGENTS_full_ledger_header_index':headings,
           'reading_scope':'named QW sources/reports, required six instruction notes, legacy action sections, current state-map/nonlinear/source boundaries; indexing does not reverify every historical programme',
           'objective_sha256':'OBJECTIVE_RECORDED_SEPARATELY_IN_COMPLETION_AUDIT',
           'source_policy':'all originals and AGENTS bytepreserved'}
    index['integration_bindings']={rel:sha(ROOT/rel) for rel in (
        'FIN_UNIFIED_EVOLUTION_SOURCE_AND_COMPOSITION_NEXT_CHAIN_20261017.md',
        'FIN_UNIFIED_EVOLUTION_AND_DH_INTEGRATION_20261017.md',
        'fin_dense_orbit_fibre_transport_review_20261017/MANIFEST.sha256',
        'fin_intermediate_clock_bridge_research_20261017/MANIFEST.sha256')}
    objective=Path('/home/krzysiek/.codex/attachments/0b4984a5-327e-48dd-aded-d2e33b8ed059/goal-objective.md')
    if objective.exists():index['objective_sha256']=sha(objective)
    indexpath=HERE/'SOURCE_INDEX.json'
    if a.check:
        old=json.loads(indexpath.read_text())
        if not objective.exists():index['objective_sha256']=old['objective_sha256']
        assert old==index
    else:indexpath.write_text(json.dumps(index,indent=2,ensure_ascii=False)+'\n')
    files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
    lines=''.join(f'{sha(p)}  {p.name}\n' for p in files)
    if a.check:assert (HERE/'MANIFEST.sha256').read_text()==lines
    else:(HERE/'MANIFEST.sha256').write_text(lines)
    print(json.dumps({'status':'PASS_UNIFIED_EVOLUTION_AUDIT_PACKAGE','entries':len(files),
                      'bounded_controls':44,'read_only':a.check,'strict_joint_origin':'NOT_ESTABLISHED'}))
if __name__=='__main__':main()
