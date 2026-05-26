from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'generated'/'p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json'
class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/'p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.py')],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2168_s1118_v1')
        self.assertTrue(d['gatekeeper_checks']['validator_exported'])
        self.assertTrue(d['gatekeeper_checks']['validator_ready'])
        self.assertFalse(d['gatekeeper_checks']['all_required_pass'])
if __name__=='__main__':
    unittest.main()
