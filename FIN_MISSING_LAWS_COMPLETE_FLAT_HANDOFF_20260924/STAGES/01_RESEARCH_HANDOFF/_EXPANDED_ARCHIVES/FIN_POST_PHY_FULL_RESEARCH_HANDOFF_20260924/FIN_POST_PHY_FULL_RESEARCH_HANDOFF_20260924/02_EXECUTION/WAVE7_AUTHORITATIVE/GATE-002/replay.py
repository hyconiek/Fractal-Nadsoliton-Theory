import json
from pathlib import Path
r=json.loads((Path(__file__).parent/'results.json').read_text())
assert r['status']=='STOP_DO_NOT_PROMOTE'
assert not r['prerequisites']['joint_action']
assert not r['prerequisites']['two_probe_universality']
print('PASS GATE-002 closed')
