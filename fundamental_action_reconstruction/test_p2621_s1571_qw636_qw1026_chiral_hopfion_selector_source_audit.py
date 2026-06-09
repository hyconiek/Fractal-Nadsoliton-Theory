from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.py"
OUT = ROOT / "generated" / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.json"
MD = ROOT / "generated" / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.md"


class P2621ChiralHopfionSelectorSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_packet_identity_and_prior_art(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2621")
        self.assertEqual(self.payload["slice_id"], "S1571")
        for item in self.payload["prior_art_files"].values():
            self.assertTrue(item["exists"])
            self.assertNotEqual(item["sha256"], "MISSING")

    def test_hopfion_energy_selector_recovers_sigma(self) -> None:
        cert = self.payload["hopfion_energy_selector_certificate"]
        self.assertTrue(cert["all_rows_recover_sigma"])
        self.assertTrue(cert["all_rows_selector_defined"])
        for row in cert["rows"]:
            self.assertEqual(row["recovered_sigma_from_delta"], row["sigma"])
            self.assertLess(row["recovery_error_abs"], 1e-12)

    def test_chiral_anomaly_selector_recovers_sigma_and_is_odd(self) -> None:
        cert = self.payload["chiral_anomaly_selector_certificate"]
        self.assertTrue(cert["all_rows_recover_sigma"])
        self.assertTrue(cert["all_rows_orientation_odd"])
        self.assertTrue(cert["all_rows_selector_defined"])

    def test_c2_torsor_has_equivariant_selector(self) -> None:
        c2 = self.payload["c2_equivariant_selector_enumeration"]
        self.assertEqual(c2["candidate_function_count"], 4)
        self.assertEqual(c2["equivariant_function_count"], 2)
        self.assertTrue(c2["canonical_selector_is_equivariant"])

    def test_p2620_not_fully_repaired_by_orientation_atom_alone(self) -> None:
        update = self.payload["p2620_conditional_update"]
        self.assertEqual(update["accepting_row_count"], 1)
        orientation_only = [
            row for row in update["rows"]
            if row["assignment"]["orientation_odd_selector_source_from_chiral_hopfion_anomaly"]
            and not row["assignment"]["nonlinear_damping_completion_source"]
        ][0]
        self.assertFalse(orientation_only["p2620_bridge_source_cut_repaired"])
        self.assertTrue(orientation_only["orientation_obstruction_repaired_conditionally"])

    def test_guard_flags_and_docs(self) -> None:
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("CONDITIONALLY_REPAIRED", self.payload["theorem_export"]["orientation_atom_status"])
        self.assertIn("P2621/S1571", MD.read_text(encoding="utf-8"))
        self.assertIn("P2621/S1571", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2621/S1571", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
