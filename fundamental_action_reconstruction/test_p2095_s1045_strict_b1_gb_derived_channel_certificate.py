from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
    ROOT / "p2095_s1045_strict_b1_gb_derived_channel_certificate.py",
]
OUT = ROOT / "generated" / "p2095_s1045_strict_b1_gb_derived_channel_certificate.json"


class TestP2095StrictB1GbDerivedChannelCertificate(unittest.TestCase):
    def test_p2095_exports_gb_derived_channel_certificate(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2095_s1045_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_B1_GB_DERIVED_CHANNEL_CERTIFICATE_WITH_TRACE__NO_FULL_4CHANNEL_INDEPENDENCE",
        )

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["gb_derived_identity_operator_level"])
        self.assertFalse(checks["full_4channel_independence_proven"])
        self.assertFalse(checks["c3_theorem_proven"])


if __name__ == "__main__":
    unittest.main()
