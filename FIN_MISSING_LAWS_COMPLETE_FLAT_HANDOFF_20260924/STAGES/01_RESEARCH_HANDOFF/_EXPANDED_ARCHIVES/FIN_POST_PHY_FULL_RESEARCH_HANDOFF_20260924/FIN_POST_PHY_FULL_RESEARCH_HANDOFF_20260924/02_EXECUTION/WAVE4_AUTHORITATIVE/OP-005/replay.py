#!/usr/bin/env python3
import json
# Example familywise budget.
deltas=[0.025,0.0125,0.0125]
assert abs(sum(deltas)-0.05)<1e-15
# Worst-case sequential first-failure probability when each reached stage uses its full conditional budget.
p_any=1.0
survive=1.0
for a in deltas:
    survive*=1-a
p_any=1-survive
assert p_any<=sum(deltas)
# Dependent marginal e-values need not multiply safely.
# B=1 with probability 1/2; e1=e2=2B.
Ee1=0.5*2+0.5*0
Ee2=Ee1
Eprod=0.5*4+0.5*0
assert Ee1==1 and Ee2==1 and Eprod==2
# Existing unmodified 5% memory + 2.5% calibration + 5% channel ledger.
unmodified_upper=0.05+0.025+0.05
assert abs(unmodified_upper-0.125)<1e-15
print(json.dumps({'allocated_deltas':deltas,'sum_bound':sum(deltas),'sequential_exact_if_full_each':p_any,'dependent_evalues':{'E1':Ee1,'E2':Ee2,'E_product':Eprod},'unmodified_subtest_union_bound':unmodified_upper},indent=2))
