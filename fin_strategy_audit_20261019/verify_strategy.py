#!/usr/bin/env python3
"""Validate the research architecture, never execute its research tasks."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs-only', action='store_true')
    args = parser.parse_args()
    ledger = json.loads((HERE/'SOURCE_LEDGER.json').read_text())
    drift = []
    for row in ledger['files']:
        p = ROOT/row['path']
        assert p.is_file(), f"Missing source: {row['path']}"
        if sha(p) != row['sha256']:
            data = p.read_bytes()
            if row['path'] == 'AGENTS.md' and len(data) >= row['bytes'] and hashlib.sha256(data[:row['bytes']]).hexdigest() == row['sha256']:
                drift.append('AGENTS has a hash-matching historical prefix; READ AND BIND THE NEW LIVE INSTRUCTIONS')
            else:
                raise AssertionError('Scientific input changed: '+row['path'])
    print('Bound source files:', len(ledger['files']), 'PASS')
    for warning in drift:
        print(warning)
    if args.inputs_only:
        return
    manifest = json.loads((ROOT/'research_tasks/TASKS.json').read_text())
    tasks = manifest['tasks']
    by_id = {t['id']: t for t in tasks}
    assert len(by_id) == len(tasks) == 5
    assert manifest['default_start'] == ['ID-001']
    assert all(t['execution_status'] == 'NOT_EXECUTED' for t in tasks)
    assert by_id['RP004-C']['dependencies'] == ['RP004-A', 'RP004-B']
    outputs = [t['output_directory'] for t in tasks]
    assert len(set(outputs)) == len(outputs)
    visiting, done = set(), set()

    def visit(key):
        assert key in by_id
        assert key not in visiting, 'Dependency cycle at '+key
        if key in done:
            return
        visiting.add(key)
        for dep in by_id[key]['dependencies']:
            visit(dep)
        visiting.remove(key)
        done.add(key)

    duplicate = json.loads((HERE/'FINAL_TASK_DUPLICATE_SEARCH.json').read_text())
    for task in tasks:
        visit(task['id'])
        p = ROOT/task['prompt']
        text = p.read_text()
        for n in range(1, 16):
            assert re.search(r'^## '+str(n)+r'\.', text, re.M), (task['id'], n)
        assert task['id'] in duplicate
        assert task['novelty_delta'] and task['stop_rule']
        for name in task['existing_inputs']:
            assert (ROOT/name).is_file(), (task['id'], name)
        for required in task['prerequisite_outputs']:
            assert task['dependencies'], (task['id'], required)
    docs = [ROOT/n for n in manifest['strategy_documents']]
    docs += list((ROOT/'research_tasks').glob('*.md'))
    docs += [HERE/'INTRINSIC_EVOLUTION_ASSESSMENT.md', HERE/'AUDIT_SCOPE.md']
    for p in docs:
        text = p.read_text()
        prose = re.sub(r'```[\s\S]*?```', '', text)
        prose = re.sub(r'`[^`\n]*`', '', prose)
        for target in re.findall(r'\]\(([^\n)]+)\)', prose):
            target = target.strip('<>').split('#')[0]
            if not target or re.match(r'https?://', target):
                continue
            assert (p.parent/target).exists(), f'Broken link: {p.name} -> {target}'
        assert not any(line.rstrip() != line for line in text.splitlines()), p
    for p in list(HERE.rglob('*.json'))+list((ROOT/'research_tasks').glob('*.json')):
        json.loads(p.read_text())
    for p in HERE.glob('*.py'):
        ast.parse(p.read_text())
    output = {
        'date': '2026-10-19', 'status': 'PASS_ARCHITECTURE_CHECKS',
        'tasks': len(tasks), 'default_start': ['ID-001'],
        'dependency_DAG': 'acyclic', 'required_sections_per_task': 15,
        'existing_inputs_and_links': 'present',
        'bound_source_files': len(ledger['files']),
        'live_context_changes': drift,
        'scientific_campaigns_executed': False,
        'future_theorems_proved_by_this_checker': False,
    }
    (HERE/'VERIFICATION.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
