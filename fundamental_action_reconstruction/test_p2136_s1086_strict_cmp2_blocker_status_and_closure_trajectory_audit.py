from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit.json"

class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/"p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit.py")],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2136_s1086_v1')
        self.assertEqual(d['result_kind'],'PASS_STRICT_CMP2_BLOCKER_STATUS_AND_CLOSURE_TRAJECTORY_AUDIT_WITH_TRACE')
        self.assertTrue(d['gatekeeper_checks']['no_toe_closure_claimed'])

if __name__=='__main__':
    unittest.main()
