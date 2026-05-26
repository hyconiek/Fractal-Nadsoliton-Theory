from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2148_s1098_strict_cmp2_external_data_delivery_blocker_report.json"
class T(unittest.TestCase):
  def test_export(self):
    subprocess.run([sys.executable,str(ROOT/"p2148_s1098_strict_cmp2_external_data_delivery_blocker_report.py")],check=True)
    d=json.loads(OUT.read_text(encoding='utf-8'))
    self.assertEqual(d['schema_version'],'p2148_s1098_v1')
    self.assertTrue(d['gatekeeper_checks']['external_data_required'])
if __name__=='__main__':
  unittest.main()
