from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.py",
    ROOT / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.py",
]
OUT = ROOT / "generated" / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.json"


class TestP2068TransportabilityImprovementLeverAudit(unittest.TestCase):
    def test_p2068_exports_ranked_levers_and_keeps_c3_open(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2068_s1018_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_TRANSPORTABILITY_IMPROVEMENT_LEVER_AUDIT_WITH_TRACE__C3_STILL_OPEN",
        )

        comps = data["transportability_components"]
        self.assertGreaterEqual(comps["transportability_gap"], 0.0)

        ranking = data["improvement_lever_ranking"]
        self.assertGreaterEqual(len(ranking), 1)
        deficits = [row["deficit_to_radius_star"] for row in ranking]
        self.assertEqual(deficits, sorted(deficits, reverse=True))

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["preconditions_ready"])
        self.assertTrue(checks["gap_nonnegative"])
        self.assertTrue(checks["ranking_nonempty"])
        self.assertTrue(checks["ranked_descending_deficit"])
        self.assertFalse(checks["c3_theorem_proven"])
        self.assertTrue(checks["no_background_globalization_claimed"])
        self.assertTrue(checks["no_tensor_component_claimed"])
        self.assertTrue(checks["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
