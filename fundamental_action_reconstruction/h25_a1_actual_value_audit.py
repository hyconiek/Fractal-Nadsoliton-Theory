#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent
out = root / 'generated'
out.mkdir(exist_ok=True)

search_files = [
    root / 'H24_A1_SOURCE_VALUE_PACKET.md',
    root / 'generated' / 'h24_a1_source_value_packet.json',
    root / 'H20_COEFFICIENT_EXPORT_SEMANTICS_PACKET.md',
    root / 'generated' / 'h20_coefficient_export_semantics_packet.json',
    root / 'generated' / 'route_a_provenance_valid_binding_instance.json',
    root / 'generated' / 'h19_first_coefficient_or_invariant_extraction.json',
    repo / 'README.md',
    repo / 'PLAN_NAPRAWA_TOE_ROBOCZY.md',
    repo / 'RAPORT_STAN_TEORII_FIN_V5_1_READINESS_2026-03-05.md',
]

hits = []
for path in search_files:
    if not path.exists():
        continue
    text = path.read_text(encoding='utf-8')
    if 'a_1' in text:
        hits.append(str(path.relative_to(repo)))

artifact = {
    'stage': 'H25',
    'target': 'a_1',
    'packet_ready_source_exists': True,
    'actual_value_witness_exists': False,
    'evaluated_numeric_value_exists': False,
    'populated_symbolic_value_exists': False,
    'partial_value_witness_exists': False,
    'files_checked_with_a1_mentions': hits,
    'frontier': 'H25_B1',
    'no_false_pass': True,
}

summary = {
    'stage': 'H25',
    'status': 'PASS_PARTIAL_A1_SOURCE_PRESENT_BUT_ACTUAL_VALUE_WITNESS_ABSENT',
    'result': 'a1_source_semantics_exist_but_no_populated_value_witness_exists',
    'frontier': [
        'H25_B1',
        'H24_B1',
        'H23_B1',
        'H22_B1',
        'H20_B1',
        'H15_B1',
        'T12_B1',
        'T2_B1',
        'C32_B2',
    ],
    'theorem_level_pass': False,
    'full_closure_pass': False,
}

(out / 'h25_a1_actual_value_audit.json').write_text(json.dumps(artifact, indent=2) + '\n', encoding='utf-8')
(out / 'h25_a1_actual_value_audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
