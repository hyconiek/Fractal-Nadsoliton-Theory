#!/usr/bin/env python3
"""Read-only repository search supporting strategic planning, not research execution."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import re
import subprocess
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
QUERIES = {
    'B01_source_action': r'same.Hessian|strict.joint.origin|nonlinear.*source|source.*nonlinear|coefficient.*map|full.action.*origin',
    'B02_kinetic_source': r'kinetic.*(source|activity)|activity.*(source|freedom)|activation.barrier|reset.*(origin|resource)|source.*(clock|generator)',
    'B03_radiation': r'RP.?004|retarded.*(flux|reference)|outward.flux|logarithmic.time|NONLINEAR_OUTWARD|flux.domination',
    'B04_stability_kernel': r'orbital.*(stability|persistence)|noncommuting|spectator.kernel|kernel.blind|active.*kernel|Phi.?6',
    'B05_Hall_IBA': r'ADJACENT|shrinking.annul|mesoscopic.*Hall|clustered.complement|unconditional.*IBA|full.orbit.Hall|shielding.weight',
    'B06_composition_geometry': r'simultaneous.*(subsystem|site|incidence)|composition.*(source|obstruction)|incidence.*source|locality.*(source|no.go)|RECORD.TO.SPATIAL',
    'B07_selector_units_bridge': r'QW.?2191|scale.charged|unit.source|role.transfer|legacy.*strict.*(completion|bridge)',
    'B08_operational_quantum': r'calibration.*(gate|source|missing)|NOT_LICENSED_BY_CURRENT_FIN_OPERATIONS|physical.*(instrument|clock)|empirical.*(absent|missing|not.supplied)',
    'B09_global_barrier': r'Gamma.?=.?B4|global.*(barrier|Target.P)|B4.*(unpaid|open|conditional)|task.?337|task.?335',
    'B10_multiscale': r'fractal.*(dynamics|no.go)|coarse.*(closure|nonclosure)|hidden.*(memory|friction)|multiscale.*(source|obstruction)',
    'B11_entropy_factor': r'Kolmogorov.?Sinai|metric entropy|Ruelle|Abramov|finite.entropy|entropy.rate.*(factor|deterministic)|Ornstein|sampling.resolution',
}


def run_query(item):
    key, pattern = item
    index = json.loads((HERE/'SEARCH_FILE_INDEX.json').read_text())['selected_paths']
    cmd = ['rg', '-l', '-i', '--threads', '1', '-e', pattern, '--']
    found = set()
    for start in range(0, len(index), 256):
        result = subprocess.run(cmd+index[start:start+256], cwd=ROOT,
                                capture_output=True, text=True, timeout=30)
        assert result.returncode in (0, 1), (key, result.stderr)
        found.update(x.removeprefix('./') for x in result.stdout.splitlines())
    files = sorted(found)
    return key, {'pattern': pattern, 'command_prefix': cmd, 'file_index': 'SEARCH_FILE_INDEX.json',
                 'batch_size': 256, 'matched_files': files, 'count': len(files),
                 'scope': 'explicit regular-file scientific-text index; see index selection; hits are not proof admission'}


def archives():
    result = subprocess.run(['rg', '--files', '-g', '*.zip'], cwd=ROOT,
                            capture_output=True, text=True, check=True)
    inventory, reads = [], []
    pattern = re.compile(r'RP.?004|retarded|outward.flux|mesoscopic|clustered.complement|nonlinear.stability', re.I)
    for name in sorted(result.stdout.splitlines()):
        path = ROOT/name
        try:
            with zipfile.ZipFile(path) as archive:
                members = archive.infolist()
                hits = [m.filename for m in members if pattern.search(m.filename)]
                record = {'path': name, 'bytes': path.stat().st_size,
                          'members': len(members), 'priority_name_hits': hits}
                inventory.append(record)
                # Read the actual four relevant RBCUT historical handoff/decision
                # texts; no extraction and no execution of archive scripts.
                if name.startswith('FIN_RBCUT_DO_WEIGHTED_DEGREE_KOMPLET_20261016/ORYGINALNE_PAKIETY/'):
                    for m in members:
                        if m.filename.endswith(('HANDOFF.md', 'NEXT_OPEN_ATOM.md', 'REPORT.md')) and m.file_size < 200000:
                            data = archive.read(m)
                            text = data.decode('utf-8', errors='replace')
                            reads.append({'archive': name, 'member': m.filename,
                                          'sha256': hashlib.sha256(data).hexdigest(),
                                          'text': text})
        except (OSError, zipfile.BadZipFile) as exc:
            inventory.append({'path': name, 'error': str(exc)})
    return inventory, reads


def main():
    results = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        pending = [pool.submit(run_query, item) for item in QUERIES.items()]
        for future in as_completed(pending):
            key, value = future.result()
            results[key] = value
            (HERE/(key+'_SEARCH.json')).write_text(json.dumps(value, indent=2)+'\n')
            print(key, value['count'], flush=True)
    (HERE/'SEARCH_RESULTS.json').write_text(json.dumps(results, indent=2)+'\n')
    inventory, reads = archives()
    (HERE/'ARCHIVE_INVENTORY.json').write_text(json.dumps(inventory, indent=2)+'\n')
    (HERE/'ARCHIVE_READS.json').write_text(json.dumps(reads, indent=2)+'\n')
    snapshot = {
        'date': '2026-10-19',
        'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'agents_sha256': hashlib.sha256((ROOT/'AGENTS.md').read_bytes()).hexdigest(),
        'query_counts': {k:v['count'] for k,v in results.items()},
        'zip_catalogues_inspected': len(inventory),
        'archive_texts_read': len(reads),
        'limits': 'No exhaustive proof replay, no proof of historical absence; not all archive payloads read. Semantic claims require separately inspected sources.'
    }
    (HERE/'SEARCH_SNAPSHOT.json').write_text(json.dumps(snapshot, indent=2)+'\n')
    print(json.dumps(snapshot, indent=2))


if __name__ == '__main__':
    main()
