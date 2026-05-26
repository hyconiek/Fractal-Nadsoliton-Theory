from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'generated'/'p2166_s1116_strict_qw2191_noncyclic_selector_source_witness.json'
class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/'p2166_s1116_strict_qw2191_noncyclic_selector_source_witness.py')],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2166_s1116_v1')
        self.assertTrue(d['gatekeeper_checks']['witness_exported'])
        self.assertTrue(d['gatekeeper_checks']['admissible_branch_selected'])
if __name__=='__main__':
    unittest.main()
