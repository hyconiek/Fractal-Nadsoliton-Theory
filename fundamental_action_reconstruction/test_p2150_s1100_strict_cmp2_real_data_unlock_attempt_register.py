from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json"

class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable, str(ROOT/"p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.py")], check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2150_s1100_v1')
        self.assertTrue(d['gatekeeper_checks']['register_exported'])

if __name__=='__main__':
    unittest.main()
