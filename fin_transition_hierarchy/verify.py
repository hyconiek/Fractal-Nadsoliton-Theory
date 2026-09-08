"""Final replay and evidence collection; never certifies prose by counting it."""
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import research
import temporal
import completion

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def canonical(value):return json.loads(json.dumps(value))


def main():
    command=[sys.executable,'-m','unittest','discover','-s','.', '-p','test_*.py','-v']
    result=subprocess.run(command,cwd=HERE,capture_output=True,text=True)
    if result.returncode:
        print(result.stderr);raise SystemExit(result.returncode)
    count=int(re.search(r'Ran (\d+) tests',result.stderr).group(1))
    replays={}
    for module,filename in [(research,'results.json'),(temporal,'temporal_results.json'),
                            (completion,'completion_results.json')]:
        expected=json.loads((HERE/filename).read_text())
        actual=canonical(module.run())
        replays[filename]=expected==actual
        if not replays[filename]:raise AssertionError(f'Archived replay differs: {filename}')
    report=(HERE/'REPORT.md').read_text()
    sections=[(int(n),int(p)) for n,p in re.findall(r'^## (\d+)\. ST(\d+)',report,re.M)]
    if sections!=[(n,8590+n) for n in range(1,31)]:
        raise AssertionError('Research-round index is incomplete or duplicated')
    for source in HERE.glob('*.py'):ast.parse(source.read_text(),filename=source.name)
    if list(HERE.glob('*.pdf')):raise AssertionError('This goal forbids PDF generation')
    source_paths=['AGENTS.md','SUMMARY_GROK.md',
        'DIAGRAMS_KERNEL_TRANSFORMATION.md','The FIN Kernel as an Unknown Mat.md',
        'FIN_Theory_Compendium_From_Fractal_Information_to_Current_Mathematics_EN.tex',
        'FIN_Post_41_50_Methodology_Correction_and_Mirror_Coupling.md',
        'FIN_Programs_31_40_Negative_Information_Coupling.md',
        'neural_coupling_investigation.py',
        'fin_replication_consistency/report.tex','fin_adaptive_30/report.tex',
        'fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md',
        'fundamental_action_reconstruction/K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md',
        'fundamental_action_reconstruction/F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md',
        'fundamental_action_reconstruction/F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md',
        'fundamental_action_reconstruction/S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md',
        'fundamental_action_reconstruction/a1_minimal_action_ansatz.py',
        'fundamental_action_reconstruction/a4_rg_emergence.py',
        'fundamental_action_reconstruction/a8_gravity_bridge.py']
    sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_paths}
    inventory=subprocess.check_output(['rg','--files'],cwd=ROOT,text=True).splitlines()
    counts={}
    for name in inventory:
        suffix=Path(name).suffix or '<none>';counts[suffix]=counts.get(suffix,0)+1
    output=dict(test_command=command,test_count=count,returncode=result.returncode,
        stdout=result.stdout,stderr=result.stderr,archived_results_replayed=replays,
        report_rounds=sections,no_pdf=True,source_sha256=sources,
        repository_inventory_count=len(inventory),inventory_by_extension=counts,
        source_scope='Repository file/state-map inventory and selected proof/source reading, not per-file historical theorem revalidation.',
        proof_scope='General theorem claims rely on the displayed mathematical proofs. Test counts and hashes are not substitutes for proofs or physical evidence.')
    (HERE/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(f'{count} scientific tests passed; all three result files reproduced exactly; 30 report sections present.')


if __name__=='__main__':main()
