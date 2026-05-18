from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1979_s929_strict_lapse_shear_provider_export_audit.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1979_s929_strict_lapse_shear_provider_export_audit.json"


class TestP1979S929(unittest.TestCase):
    def test_lapse_shear_provider_export_audit(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1979")
        self.assertEqual(payload["stage_id"], "S929")
        self.assertEqual(payload["result_kind"], "PASS_CURRENT_EXPORT_NONAVAILABILITY_AUDIT")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["depends_on"]["p1978_bounded_obstruction_present"])
        self.assertEqual(payload["required_non_energy_neutral_provider_signature"]["required_U00"], "-sigma1**2 - sigma1*sigma2 - sigma2**2")
        self.assertEqual(payload["required_non_energy_neutral_provider_signature"]["required_spatial_trace"], "-3*(sigma1**2 + sigma1*sigma2 + sigma2**2)")
        self.assertTrue(payload["gatekeeper_checks"]["no_exported_lapse_shear_provider_certificate_in_current_registries"])
        self.assertTrue(payload["gatekeeper_checks"]["bounded_current_export_audit_passed"])
        self.assertIn("not a theorem", payload["theorem_export"]["not_a_no_go_clause"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
