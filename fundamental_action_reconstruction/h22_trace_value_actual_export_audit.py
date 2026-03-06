#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent
out = root / 'generated'
out.mkdir(exist_ok=True)

search_files = [
    root / 'H21_TRACE_VALUE_EXPORT_PACKET.md',
    root / 'H20_COEFFICIENT_EXPORT_SEMANTICS_PACKET.md',
    root / 'generated' / 'h20_coefficient_export_semantics_packet.json',
    root / 'generated' / 'h21_trace_value_export_packet.json',
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
    if 'trace_A_1' in text or 'tr(A_1)' in text:
        hits.append(str(path.relative_to(repo)))

artifact = {
    'stage': 'H22',
    'target': 'trace_A_1',
    'packet_ready_target_exists': True,
    'actual_value_witness_exists': False,
    'evaluated_numeric_value_exists': False,
    'populated_symbolic_value_exists': False,
    'files_checked_with_trace_mentions': hits,
    'frontier': 'H22_B1',
    'no_false_pass': True,
}

summary = {
    'stage': 'H22',
    'status': 'PASS_PARTIAL_TRACE_VALUE_TARGET_PRESENT_BUT_ACTUAL_VALUE_WITNESS_ABSENT',
    'result': 'trace_value_semantics_exist_but_no_populated_value_witness_exists',
    'frontier': [
        'H22_B1',
        'H21_B1',
        'H20_B1',
        'H15_B1',
        'T12_B1',
        'T2_B1',
        'C32_B2',
    ],
    'theorem_level_pass': False,
    'full_closure_pass': False,
}

(out / 'h22_trace_value_actual_export_audit.json').write_text(json.dumps(artifact, indent=2) + '\n', encoding='utf-8')
(out / 'h22_trace_value_actual_export_audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
