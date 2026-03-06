from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
cert = json.loads((root / 'generated' / 'axiom_lane_mode_scaffold_compatibility_certificate.json').read_text(encoding='ascii'))
summary = json.loads((root / 'generated' / 'ax5_axiom_lane_mode_scaffold_compatibility_audit_summary.json').read_text(encoding='ascii'))
print(root / 'generated' / 'ax5_axiom_lane_mode_scaffold_compatibility_audit_summary.json')
