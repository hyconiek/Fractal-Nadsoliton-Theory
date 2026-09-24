#!/usr/bin/env python3
import json
x={"1.5": {"degeneracy_counts": {"1": 366, "2": 126, "4": 3, "other": 0}, "median_gap": 0.14058028044428283, "median_n95": 24.0, "max_n95_capped_1001": 1001}, "2.0": {"degeneracy_counts": {"1": 366, "2": 126, "4": 3, "other": 0}, "median_gap": 0.33739267306628307, "median_n95": 10.0, "max_n95_capped_1001": 130}, "3.0": {"degeneracy_counts": {"1": 366, "2": 126, "4": 3, "other": 0}, "median_gap": 0.8997137948433647, "median_n95": 4.0, "max_n95_capped_1001": 49}}
print(json.dumps(x,indent=2,sort_keys=True))
assert x['2.0']['degeneracy_counts']['1']==366
assert x['2.0']['degeneracy_counts']['2']==126
assert x['2.0']['degeneracy_counts']['4']==3
assert x['2.0']['median_n95']==10.0
