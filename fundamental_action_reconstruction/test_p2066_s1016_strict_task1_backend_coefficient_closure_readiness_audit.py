import json
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SCRIPT = ROOT / "p2066_s1016_strict_task1_backend_coefficient_closure_readiness_audit.py"
OUT = GEN / "p2066_s1016_strict_task1_backend_coefficient_closure_readiness_audit.json"


class TestP2066(unittest.TestCase):
    def test_packet(self):
        subprocess.run(["python3", str(SCRIPT)], check=True)
        self.assertTrue(OUT.exists())
        data = json.loads(OUT.read_text(encoding="utf-8"))

        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_TASK1_BACKEND_COEFFICIENT_CLOSURE_READINESS_AUDIT__QUOTIENT_ONLY_STILL_OPEN",
        )
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")

        grep_block = data["content_grep_evidence"]
        self.assertGreaterEqual(grep_block["match_count"], 1)

        readiness = data["closure_readiness"]
        self.assertTrue(readiness["backend_coefficients_present"])
        self.assertFalse(readiness["independent_a_GB_identified"])
        self.assertFalse(readiness["strict_task1_full_closure_ready"])

        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["no_false_pass"])
        self.assertTrue(checks["c3_transport_theorem_open"])


if __name__ == "__main__":
    unittest.main()
