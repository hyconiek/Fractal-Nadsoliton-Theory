from __future__ import annotations
import json, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
G=ROOT/'generated'
class T(unittest.TestCase):
    def test_packets(self):
        subprocess.run([sys.executable,str(ROOT/'p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.py')],check=True)
        subprocess.run([sys.executable,str(ROOT/'p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.py')],check=True)
        d1=json.loads((G/'p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.json').read_text())
        d2=json.loads((G/'p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.json').read_text())
        self.assertEqual(d1['schema_version'],'p2154_s1104_v1')
        self.assertEqual(d2['schema_version'],'p2155_s1105_v1')
        self.assertTrue(d1['gatekeeper_checks']['milestone_packet_exported'])
        self.assertTrue(d2['gatekeeper_checks']['gap_packet_exported'])
if __name__=='__main__':
    unittest.main()
