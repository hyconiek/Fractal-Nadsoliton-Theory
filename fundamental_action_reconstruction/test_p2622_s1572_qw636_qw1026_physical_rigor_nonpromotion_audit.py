from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.py"
OUT = ROOT / "generated" / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.json"
MD = ROOT / "generated" / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.md"


class P2622QW636QW1026PhysicalRigorNonpromotionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_packet_identity_and_prior_art_files(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2622")
        self.assertEqual(self.payload["slice_id"], "S1572")
        for item in self.payload["prior_art_files"].values():
            self.assertTrue(item["exists"])
            self.assertNotEqual(item["sha256"], "MISSING")

    def test_q636_is_covariant_and_open_chain_phase_is_gauge(self) -> None:
        cert = self.payload["q636_covariance_and_gauge_certificate"]
        self.assertLess(cert["parity_covariance_max_defect"], 1e-12)
        self.assertLess(cert["sorted_spectrum_sigma_plus_vs_minus_max_defect"], 1e-12)
        self.assertLess(cert["open_chain_gauge_equivalence_max_defect"], 1e-12)
        self.assertLess(cert["open_chain_eigenvalue_max_defect"], 1e-12)
        self.assertIn("CONDITIONAL_DIAGNOSTIC", cert["methodological_verdict"])

    def test_q1026_anomaly_sign_flips_with_unsourced_gamma5_origin(self) -> None:
        cert = self.payload["q1026_gamma5_convention_certificate"]
        self.assertAlmostEqual(cert["opposite_gamma5_sign_sum_real"], 0.0, places=9)
        self.assertAlmostEqual(cert["opposite_gamma5_sign_sum_imag"], 0.0, places=9)
        self.assertTrue(cert["shifted_gamma5_equals_minus_original"])
        self.assertTrue(cert["sign_flips_under_allowed_alternating_origin_change"])
        self.assertIn("CONVENTION_DEPENDENT", cert["methodological_verdict"])

    def test_prior_art_alone_is_not_accepted_as_selector_source(self) -> None:
        self.assertFalse(self.payload["accepted_orientation_source_from_qw636_qw1026_prior_art_alone"])
        rows = self.payload["strict_acceptance_table"]
        self.assertTrue(rows)
        self.assertTrue(all(not row["accepted_as_orientation_odd_selector_source_now"] for row in rows))
        self.assertEqual(self.payload["p2620_prior_art_update"]["accepting_row_count"], 0)

    def test_guard_flags_and_docs(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2622/S1572", MD.read_text(encoding="utf-8"))
        self.assertIn("P2622/S1572", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2622/S1572", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
