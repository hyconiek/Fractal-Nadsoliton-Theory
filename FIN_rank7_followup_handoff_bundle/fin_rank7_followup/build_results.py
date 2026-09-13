#!/usr/bin/env python3
# Explicit regeneration entry point. It never rewrites imported inputs.
from pathlib import Path
import json, datetime
ROOT=Path(__file__).resolve().parent
out=ROOT/'results'/'build_snapshot.json'
obj={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
     'task_counts':{},'claim_counts':{}}
for t in json.loads((ROOT/'TASKS.json').read_text()): obj['task_counts'][t['execution_status']]=obj['task_counts'].get(t['execution_status'],0)+1
for c in json.loads((ROOT/'CLAIMS.json').read_text()): obj['claim_counts'][c['status']]=obj['claim_counts'].get(c['status'],0)+1
out.write_text(json.dumps(obj,indent=2)+'\n')
print(out)
