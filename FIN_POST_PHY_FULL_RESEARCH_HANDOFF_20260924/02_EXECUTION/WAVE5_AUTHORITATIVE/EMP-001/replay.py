#!/usr/bin/env python3
import json
required=['preparation_reset_fidelity','independent_detector_memory_bound','complete_outcome_map','control_map_calibration_eta_and_t1','efficiency_confusion_model','calibration_drift_bound','per_attempt_reset_costs','raw_data_custody_and_split_provenance']
platform={}
missing=[k for k in required if k not in platform]
assert missing==required
costs={'carrier_specific_control':31,'carrier_one_state_top2':435,'carrier_full_per_state':612,'channel05':31299102,'op002':471000000,'memory_upper':4294967296}
assert costs['carrier_specific_control']<costs['carrier_one_state_top2']<costs['channel05']<costs['op002']<costs['memory_upper']
print(json.dumps({'status':'NOT_READY','missing':missing,'cost_order':costs},indent=2))
