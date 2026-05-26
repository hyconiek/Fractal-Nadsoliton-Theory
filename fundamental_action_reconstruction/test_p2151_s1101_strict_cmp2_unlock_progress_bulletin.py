from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2151_s1101_strict_cmp2_unlock_progress_bulletin.json"
class T(unittest.TestCase):
  def test_export(self):
    subprocess.run([sys.executable,str(ROOT/"p2151_s1101_strict_cmp2_unlock_progress_bulletin.py")],check=True)
    d=json.loads(OUT.read_text(encoding='utf-8'))
    self.assertEqual(d['schema_version'],'p2151_s1101_v1')
    self.assertTrue(d['gatekeeper_checks']['bulletin_exported'])
if __name__=='__main__':
  unittest.main()
