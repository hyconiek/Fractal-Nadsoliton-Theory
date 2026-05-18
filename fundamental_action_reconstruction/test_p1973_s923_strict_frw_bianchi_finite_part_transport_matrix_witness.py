from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1973_s923_strict_frw_bianchi_finite_part_transport_matrix_witness.json"


class TestP1973S923(unittest.TestCase):
    def test_local_transport_matrix_witness(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1973")
        self.assertEqual(payload["stage_id"], "S923")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(payload["local_result_kind"], "PASS_ZERO_LOCAL_TRANSPORT_RESIDUAL")
        self.assertTrue(payload["symbolic_transport"]["transport_residual_is_zero"])
        self.assertTrue(payload["symbolic_transport"]["matrix_residual_zero"])
        self.assertTrue(payload["renormalization_lock_from_p1972"]["zero_vector_transport_preserved"])
        self.assertTrue(payload["numeric_transport_replay"]["all_scipy_expm_rows_match_closed_form"])
        self.assertIn("not global background-independence closure", payload["false_pass_guard"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
