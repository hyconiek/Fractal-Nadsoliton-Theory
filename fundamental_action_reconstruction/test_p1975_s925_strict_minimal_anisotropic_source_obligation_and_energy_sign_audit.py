from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json"


class TestP1975S925(unittest.TestCase):
    def test_minimal_source_obligation_and_sign_audit(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1975")
        self.assertEqual(payload["stage_id"], "S925")
        self.assertEqual(payload["result_kind"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["symbolic_checks"]["cancellation_if_source_admitted"])
        self.assertEqual(payload["symbolic_checks"]["post_source_residual"], ["0", "0", "0", "0"])
        self.assertTrue(payload["symbolic_checks"]["source_frw_limit_zero"])
        self.assertTrue(payload["symbolic_checks"]["rho_required_equals_minus_q_shear"])
        self.assertTrue(payload["symbolic_checks"]["q_shear_positive_definite"])
        self.assertTrue(payload["gatekeeper_checks"]["required_energy_density_negative_for_nonzero_shear"])
        self.assertFalse(payload["gatekeeper_checks"]["strict_source_derivation_exported"])
        self.assertIn("REQUIRED_RHO_EQUALS_NEGATIVE_SHEAR_QUADRATIC", payload["obstruction_tags"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
