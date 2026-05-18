from __future__ import annotations
import json, subprocess, unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parent/"p1991_s941_strict_augmented_provider_channel_matrix_witness.py"
OUT=Path(__file__).resolve().parent/"generated"/"p1991_s941_strict_augmented_provider_channel_matrix_witness.json"

class TestP1991S941(unittest.TestCase):
    def test_channel_matrix(self)->None:
        subprocess.run(["python3",str(SCRIPT)],check=True)
        p=json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(p["packet_id"],"P1991")
        self.assertEqual(p["result_kind"],"PASS_CHANNEL_MATRIX_OBSTRUCTION_WITNESS")
        self.assertTrue(p["gatekeeper_checks"]["matrix_rank_aug_exceeds_rankA_or_not_full_solve"])
        self.assertTrue(p["gatekeeper_checks"]["least_squares_residual_positive"])

if __name__=="__main__":
    unittest.main()
