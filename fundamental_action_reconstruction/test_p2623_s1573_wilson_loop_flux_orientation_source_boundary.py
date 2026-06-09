from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.py"
OUT = ROOT / "generated" / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.json"
MD = ROOT / "generated" / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.md"


class P2623WilsonLoopFluxOrientationSourceBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_packet_identity_and_content_first_grep(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2623")
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertTrue(all(data["count"] > 0 for data in audit["patterns"].values()))

    def test_wilson_loop_is_gauge_invariant_and_reversal_odd(self) -> None:
        cert = self.payload["wilson_flux_certificate"]
        self.assertLess(cert["gauge_invariance_defect_abs"], 1e-12)
        self.assertLess(cert["reversal_conjugacy_defect_abs"], 1e-12)
        self.assertTrue(cert["orientation_sign_flips"])
        self.assertTrue(cert["selector_from_oriented_cycle_defined"])

    def test_unoriented_wilson_data_are_sign_blind(self) -> None:
        cert = self.payload["wilson_flux_certificate"]
        self.assertTrue(cert["unoriented_class_contains_both_signs"])
        self.assertFalse(cert["selector_from_unoriented_cycle_defined"])

    def test_acceptance_lattice_requires_all_three_sources(self) -> None:
        lattice = self.payload["orientation_source_acceptance_lattice"]
        self.assertEqual(len(lattice["rows"]), 8)
        self.assertEqual(lattice["accepting_row_count"], 1)
        self.assertFalse(lattice["current_orientation_repaired"])
        accepting = [row for row in lattice["rows"] if row["orientation_odd_selector_source_repaired"]]
        self.assertEqual(len(accepting), 1)
        self.assertTrue(all(accepting[0]["assignment"].values()))

    def test_p2620_and_ltotal_remain_closed(self) -> None:
        self.assertEqual(self.payload["p2620_update"]["accepting_row_count"], 0)
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2623/S1573", MD.read_text(encoding="utf-8"))
        self.assertIn("P2623/S1573", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2623/S1573", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
