from __future__ import annotations
import json, subprocess, tempfile, unittest
from pathlib import Path
SCRIPT = Path(__file__).resolve().parent / 'p1311_qw2191_r27_nb1_archival_and_external_disclosure_packet_checkpoint.py'

class TestP1311(unittest.TestCase):
    def test_happy_path(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p1310 = d/'p1310.json'; out=d/'p1311.json'; rel=d/'R8.md'
            rel.write_text('ok', encoding='utf-8')
            p1310.write_text(json.dumps({'next_priority':'R27_NB1_ARCHIVAL_AND_EXTERNAL_DISCLOSURE_PACKET','r26_independent_replay_audit':{'status':'AUDIT_COMPLETE'}}), encoding='utf-8')
            subprocess.run(['python3', str(SCRIPT), '--p1310', str(p1310), '--release', str(rel), '--out', str(out)], check=True)
            payload=json.loads(out.read_text(encoding='utf-8'))
            self.assertEqual(payload['disclosure']['status'], 'READY')

    def test_requires_audit_complete(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            p1310 = d/'p1310.json'; out=d/'p1311.json'; rel=d/'R8.md'
            rel.write_text('ok', encoding='utf-8')
            p1310.write_text(json.dumps({'next_priority':'R27_NB1_ARCHIVAL_AND_EXTERNAL_DISCLOSURE_PACKET','r26_independent_replay_audit':{'status':'BLOCKED'}}), encoding='utf-8')
            p = subprocess.run(['python3', str(SCRIPT), '--p1310', str(p1310), '--release', str(rel), '--out', str(out)], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn('AUDIT_COMPLETE', p.stderr)

if __name__ == '__main__':
    unittest.main()
