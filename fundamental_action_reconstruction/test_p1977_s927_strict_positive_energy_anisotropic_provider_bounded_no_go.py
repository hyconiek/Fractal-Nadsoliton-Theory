from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1977_s927_strict_positive_energy_anisotropic_provider_bounded_no_go.json"


class TestP1977S927(unittest.TestCase):
    def test_positive_energy_bounded_no_go(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1977")
        self.assertEqual(payload["stage_id"], "S927")
        self.assertEqual(payload["result_kind"], "PASS_BOUNDED_NO_GO")
        self.assertEqual(payload["status"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["symbolic_core"]["rho_required_equals_minus_Q_shear"])
        self.assertTrue(payload["symbolic_core"]["Q_positive_definite"])
        self.assertTrue(payload["gatekeeper_checks"]["p1976_current_basis_has_no_exported_provider"])
        self.assertTrue(payload["gatekeeper_checks"]["all_nonzero_shear_samples_violate_positive_energy"])
        self.assertTrue(payload["gatekeeper_checks"]["bounded_no_go_passed"])
        self.assertIn("POSITIVE_ENERGY_MINIMAL_ANISOTROPIC_PROVIDER_NO_GO_UNDER_CURRENT_BASIS", payload["obstruction_tags"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
