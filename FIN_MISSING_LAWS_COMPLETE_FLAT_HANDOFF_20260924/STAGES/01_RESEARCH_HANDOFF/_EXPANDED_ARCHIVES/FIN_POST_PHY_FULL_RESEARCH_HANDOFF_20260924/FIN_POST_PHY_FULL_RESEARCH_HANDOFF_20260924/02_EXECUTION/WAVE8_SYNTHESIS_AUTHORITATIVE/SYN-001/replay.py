import json
from pathlib import Path
r=json.loads((Path(__file__).parent/'results.json').read_text())
assert r['global_pass'] is False
assert r['status'].startswith('PARTIAL_REVIEW_READY')
assert len(r['missing_for_full_audit'])==3
print('PASS SYN-001 fail-closed provenance audit')
