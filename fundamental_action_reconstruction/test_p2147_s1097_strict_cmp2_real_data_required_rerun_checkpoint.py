from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json"
class T(unittest.TestCase):
  def test_export(self):
    subprocess.run([sys.executable,str(ROOT/"p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.py")],check=True)
    d=json.loads(OUT.read_text(encoding='utf-8'))
    self.assertEqual(d['schema_version'],'p2147_s1097_v1')
    self.assertIn(d['result_kind'],{'PASS_STRICT_CMP2_REAL_DATA_REQUIRED_RERUN_CHECKPOINT_WITH_TRACE','OPEN_STRICT_CMP2_REAL_DATA_REQUIRED_RERUN_CHECKPOINT_BLOCKED'})
    self.assertTrue(d['gatekeeper_checks']['no_toe_closure_claimed'])
if __name__=='__main__':
  unittest.main()
