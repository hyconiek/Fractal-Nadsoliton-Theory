from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2146_s1096_strict_cmp2_commitment_execution_audit.json"

class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/"p2146_s1096_strict_cmp2_commitment_execution_audit.py")],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2146_s1096_v1')
        self.assertEqual(d['result_kind'],'PASS_STRICT_CMP2_COMMITMENT_EXECUTION_AUDIT_WITH_TRACE')
        self.assertTrue(d['gatekeeper_checks']['audit_exported'])

if __name__=='__main__':
    unittest.main()
