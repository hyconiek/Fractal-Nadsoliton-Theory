#!/usr/bin/env python3
import json
prereq={'SRC-003':'classical orthogonal gauge fixing','OP-004':'classical detector-memory alias'}
nonclassical_required=['noncommuting_effect_algebra','sourced_composition_rule','operational_independence']
present=[]
assert not present
print(json.dumps({'prerequisites':prereq,'required_new_objects':nonclassical_required,'present':present,'decision':'STOP_DO_NOT_PROMOTE'},indent=2))
