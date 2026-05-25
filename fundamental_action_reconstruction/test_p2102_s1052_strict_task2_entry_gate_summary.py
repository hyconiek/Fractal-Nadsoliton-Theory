from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    ROOT / "p1950_s900_strict_renormalization_exact_integration.py",
    ROOT / "p2094_s1044_strict_b1_quotient_renormalization_rank_repair.py",
    ROOT / "p2095_s1045_strict_b1_gb_derived_channel_certificate.py",
    ROOT / "p2096_s1046_strict_b1_renormalization_closure_contract.py",
    ROOT / "p2097_s1047_strict_b1_quotient_closure_stability_minigrid.py",
    ROOT / "p2098_s1048_strict_precutkosky_readiness_contract.py",
    ROOT / "p2099_s1049_strict_u1_same_scheme_lock_witness.py",
    ROOT / "p2100_s1050_strict_u2_phase_space_quadrature_witness.py",
    ROOT / "p2101_s1051_strict_u3_residue_positivity_uncertainty_witness.py",
    ROOT / "p2102_s1052_strict_task2_entry_gate_summary.py",
]
OUT = ROOT / "generated" / "p2102_s1052_strict_task2_entry_gate_summary.json"


class TestP2102StrictTask2EntryGateSummary(unittest.TestCase):
    def test_p2102_exports_gate_summary(self) -> None:
        for script in SCRIPTS:
            subprocess.run([sys.executable, str(script)], check=True)

        data = json.loads(OUT.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], "p2102_s1052_v1")
        self.assertEqual(data["status"], "OPEN_PARTIAL_PROGRESS_WITH_TRACE")
        self.assertEqual(
            data["result_kind"],
            "PASS_STRICT_TASK2_ENTRY_GATE_SUMMARY_WITH_TRACE__U1_U2_U3_READY",
        )
        checks = data["gatekeeper_checks"]
        self.assertTrue(checks["u1_ready"])
        self.assertTrue(checks["u2_ready"])
        self.assertTrue(checks["u3_ready"])
        self.assertTrue(checks["gate_ready"])
        self.assertFalse(checks["full_cutkosky_closure_proven"])


if __name__ == "__main__":
    unittest.main()
