import json
from pathlib import Path
r=json.loads((Path(__file__).parent/'results.json').read_text())
assert r['automatic_downstream_campaign'] is False
assert len(r['recommendations'])<=3
assert all(x['priority_score']>0 for x in r['recommendations'])
print('PASS SYN-002 handoff decision')
