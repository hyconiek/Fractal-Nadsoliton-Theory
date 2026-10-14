from __future__ import annotations
import json
from pathlib import Path
from theory_core import fold_generator_metrics

LAMBDAS = [0.0, 0.5, 1.0, -0.5]
NS = [120, 240, 480, 960, 1920, 3840]

rows = []
for lam in LAMBDAS:
    for N in NS:
        r = fold_generator_metrics(N, lam)
        if 'lambda_' in r:
            r['lambda'] = r.pop('lambda_')
        rows.append(r)

out = Path(__file__).resolve().parent / 'TH-005' / 'FOLD_GENERATOR_DIAGNOSTICS.json'
out.write_text(json.dumps(rows, indent=2, sort_keys=True) + '\n')
print(json.dumps({'rows': len(rows), 'output': str(out)}, indent=2))
