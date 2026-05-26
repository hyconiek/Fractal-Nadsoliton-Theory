from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"generated"/"p2145_s1095_strict_cmp2_next_action_commitment_packet.json"

class T(unittest.TestCase):
    def test_export(self):
        subprocess.run([sys.executable,str(ROOT/"p2145_s1095_strict_cmp2_next_action_commitment_packet.py")],check=True)
        d=json.loads(OUT.read_text(encoding='utf-8'))
        self.assertEqual(d['schema_version'],'p2145_s1095_v1')
        self.assertEqual(d['result_kind'],'PASS_STRICT_CMP2_NEXT_ACTION_COMMITMENT_PACKET_WITH_TRACE')
        self.assertTrue(d['gatekeeper_checks']['commitment_exported'])

if __name__=='__main__':
    unittest.main()
