from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SCRIPTS=[
 ROOT/"p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.py",
 ROOT/"p2127_s1077_strict_cmp2_bootstrap_backend_evidence_stresstest.py",
 ROOT/"p2131_s1081_strict_cmp2_support_expansion_replay_audit.py",
]
OUT=ROOT/"generated"/"p2131_s1081_strict_cmp2_support_expansion_replay_audit.json"

class TestP2131(unittest.TestCase):
    def test_export(self):
        for s in SCRIPTS:
            subprocess.run([sys.executable,str(s)],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2131_s1081_v1')
        self.assertEqual(d['result_kind'],'PASS_STRICT_CMP2_SUPPORT_EXPANSION_REPLAY_AUDIT_WITH_TRACE')
        self.assertTrue(d['gatekeeper_checks']['replay_executed'])
        self.assertIn('x2',d['support_expansion_replay_audit']['results'])

if __name__=='__main__':
    unittest.main()
