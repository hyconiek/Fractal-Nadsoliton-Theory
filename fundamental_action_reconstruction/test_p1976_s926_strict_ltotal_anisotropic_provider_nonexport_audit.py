from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.py"
OUT = Path(__file__).resolve().parent / "generated" / "p1976_s926_strict_ltotal_anisotropic_provider_nonexport_audit.json"


class TestP1976S926(unittest.TestCase):
    def test_ltotal_anisotropic_provider_nonexport_audit(self) -> None:
        subprocess.run(["python3", str(SCRIPT)], check=True)
        payload = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(payload["packet_id"], "P1976")
        self.assertEqual(payload["stage_id"], "S926")
        self.assertEqual(payload["result_kind"], "OPEN_OBSTRUCTION_WITH_TRACE")
        self.assertTrue(payload["target_obligation_from_p1975"]["rho_required_equals_minus_q_shear"])
        self.assertGreater(payload["audit_scope"]["terms_audited"], 0)
        self.assertEqual(payload["feature_replay"]["candidate_count"], 0)
        self.assertEqual(payload["feature_replay"]["total_shear_token_hits"], 0)
        self.assertTrue(payload["gatekeeper_checks"]["no_explicit_anisotropic_provider_in_current_registries"])
        self.assertTrue(payload["gatekeeper_checks"]["nonexport_not_universal_no_go"])
        self.assertIn("STRICT_LTOTAL_ANISOTROPIC_PROVIDER_NOT_EXPORTED", payload["obstruction_tags"])
        self.assertIn("ToE closure", payload["theorem_export"]["not_licensed"])


if __name__ == "__main__":
    unittest.main()
