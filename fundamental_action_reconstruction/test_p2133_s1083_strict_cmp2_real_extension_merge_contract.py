from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2133_s1083_strict_cmp2_real_extension_merge_contract.json"

class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/"p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py")],check=True)
        subprocess.run([sys.executable,str(ROOT/"p2133_s1083_strict_cmp2_real_extension_merge_contract.py")],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2133_s1083_v1')
        self.assertEqual(d['status'],'OPEN_PARTIAL_PROGRESS_WITH_TRACE')
        self.assertIn(d['result_kind'],{'PASS_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_WITH_TRACE','OPEN_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_BLOCKED'})

if __name__=='__main__':
    unittest.main()
