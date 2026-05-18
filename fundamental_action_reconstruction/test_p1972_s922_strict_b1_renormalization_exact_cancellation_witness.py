from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1972_s922_strict_b1_renormalization_exact_cancellation_witness.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1972_s922_strict_b1_renormalization_exact_cancellation_witness.json"


class TestP1972S922(unittest.TestCase):
    def test_exact_b1_counterterm_cancellation_witness(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        required_schema = {
            "ledger_id",
            "produced_by",
            "background_family_id",
            "index_convention_id",
            "boundary_clause_id",
            "component_basis",
            "result_kind",
            "residual_vector",
            "obstruction_tags",
            "timestamp_utc",
        }
        self.assertTrue(required_schema.issubset(payload))
        self.assertEqual(payload["packet_id"], "P1972")
        self.assertEqual(payload["stage_id"], "S922")
        self.assertEqual(payload["result_kind"], "PASS_ZERO")
        self.assertEqual(payload["residual_vector"], ["0", "0", "0", "0"])
        self.assertEqual(payload["obstruction_tags"], [])
        self.assertTrue(payload["algebraic_cancellation"]["residual_is_zero"])
        self.assertTrue(payload["gatekeeper_checks"]["p1853_consistency_all_match"])
        self.assertTrue(payload["gatekeeper_checks"]["rational_probe_residuals_zero"])
        self.assertIn("only to exact B1", payload["false_pass_guard"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
