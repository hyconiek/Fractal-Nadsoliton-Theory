#!/usr/bin/env python3
import json
checks={'noncommuting_instruments':False,'sourced_composite_rule':False,'classical_memory_excluded':False,'affine_preparation_novelty':False}
assert not all(checks.values())
print(json.dumps({'checks':checks,'decision':'STOP_DO_NOT_PROMOTE'},indent=2))
