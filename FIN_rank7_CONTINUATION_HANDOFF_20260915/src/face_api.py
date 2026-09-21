"""R7P-024 read-only interface for the actually certified face statements."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]

def verify_face_certificate():
 p=ROOT/'certificates/R7P-017_021_face_certificate.json'
 q=ROOT/'certificates/R7P-023_parity_first_order.json'
 if not p.exists() or not q.exists():return {'status':'MISSING'}
 a=json.loads(p.read_text());b=json.loads(q.read_text())
 ok=(a['face_resolvent']['all_positive'] and a['extreme_face']['all_positive']
     and a['unique_minimum']['ok'] and b['both_strictly_negative'])
 return {'status':'PASS' if ok else 'FAIL',
   'certified':{
     'corrected_resolvent_positive_on_r_0_1':a['face_resolvent']['all_positive'],
     'unique_resolvent_minimum_in_box':a['unique_minimum']['root_box'] if a['unique_minimum']['ok'] else None,
     'extreme_face_second_curvature_le_sigma':a['extreme_face']['all_positive'],
     'finite_equality_absent':not a['equality_locus']['finite_equality'],
     'parity_first_order_shifts_negative':b['both_strictly_negative']},
   'not_certified':['global boundary-Ising curvature','intraparity W closure','off-face s4/s5 ceiling','reoptimized parity-envelope slope']}
if __name__=='__main__':print(json.dumps(verify_face_certificate(),indent=2))
