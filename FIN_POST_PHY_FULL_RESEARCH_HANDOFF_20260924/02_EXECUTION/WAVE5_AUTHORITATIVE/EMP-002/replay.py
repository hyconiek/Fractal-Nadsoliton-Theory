#!/usr/bin/env python3
import math,json
C=0.07654211034
# Equal-prior Chernoff upper bound 0.5 exp(-n C) <= .05 => n >= log(10)/C
n=math.ceil(math.log(10)/C)
assert n==31
counts={'specific_control':n,'one_state_top2':435,'full_12_state':612*12,'channel05':31299102,'op002':471000000,'memory_upper':4294967296}
assert counts['specific_control']<counts['one_state_top2']<counts['full_12_state']<counts['channel05']<counts['op002']<counts['memory_upper']
print(json.dumps({'chernoff_required_exits':n,'ordered_cost_scales':counts},indent=2))
