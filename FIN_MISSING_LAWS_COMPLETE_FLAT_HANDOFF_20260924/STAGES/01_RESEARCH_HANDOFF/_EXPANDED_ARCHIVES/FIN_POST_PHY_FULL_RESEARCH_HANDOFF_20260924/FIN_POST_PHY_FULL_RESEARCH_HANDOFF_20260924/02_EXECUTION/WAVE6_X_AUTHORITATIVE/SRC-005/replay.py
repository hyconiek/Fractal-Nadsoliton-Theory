from pathlib import Path
import json
r=json.loads((Path(__file__).parent/'results.json').read_text())
assert r['status']=='NO_ADMISSIBLE_SOURCE'
assert r['candidate_count_admitted']==0
assert r['repo_head'].startswith('d4c4a0ac')
print('PASS SRC-005 admission record')
