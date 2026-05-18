from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json"


class TestP1974S924(unittest.TestCase):
    def test_anisotropic_eom_obstruction_witness(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1974")
        self.assertEqual(payload["stage_id"], "S924")
        self.assertEqual(payload["result_kind"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(all(payload["polynomial_nonzero_flags"]))
        self.assertEqual(payload["isotropic_limit_residual_vector"], ["0", "0", "0", "0"])
        self.assertTrue(payload["isotropic_limit_zero"])
        self.assertTrue(payload["gatekeeper_checks"]["generic_anisotropic_residual_nonzero"])
        self.assertTrue(payload["gatekeeper_checks"]["numeric_samples_nonzero"])
        self.assertIn("SCALAR_TRANSPORT_INSUFFICIENT_FOR_TENSOR_EOM", payload["obstruction_tags"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
