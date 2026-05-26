from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json"


class TestP2144StrictCmp2RealCiStabilityInterpretationGate(unittest.TestCase):
    def test_p2144_exports_interpretation_gate(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.py")], check=True)
        d = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(d["schema_version"], "p2144_s1094_v1")
        self.assertIn(d["result_kind"], {
            "PASS_STRICT_CMP2_REAL_CI_STABILITY_INTERPRETATION_GATE_WITH_TRACE",
            "OPEN_STRICT_CMP2_REAL_CI_STABILITY_INTERPRETATION_GATE_BLOCKED",
        })
        self.assertTrue(d["gatekeeper_checks"]["no_toe_closure_claimed"])


if __name__ == "__main__":
    unittest.main()
