from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.py"
OUT = ROOT / "generated" / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.json"
MD = ROOT / "generated" / "p2664_s1614_boundary_phase_sector_selector_variational_no_go_audit.md"


class P2664BoundaryPhaseSectorSelectorVariationalNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_and_p2663_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in ("sector_selector_content", "variational_phase_content", "entropy_bit_target_content", "nonclosure_guard_content"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2663_nonexact_bit_carrier_found"])
        self.assertTrue(upstream["p2663_unique_target_not_derived"])
        self.assertTrue(upstream["p2663_boundary_phase_bit_target_not_exported"])
        self.assertTrue(upstream["p2663_qw2191_not_discharged"])

    def test_variational_witness_blocks_local_selector(self) -> None:
        witness = self.payload["sector_selector_witness"]
        self.assertEqual(witness["triangle_flat_cochain_count"], 32)
        self.assertEqual(witness["square_holonomy_sector_counts_after_flatness"], {"0": 16, "1": 16})
        self.assertEqual(witness["exact_coboundary_sector_counts"], {"0": 16, "1": 0})
        self.assertFalse(witness["local_even_class_selects_nonexact_sector"])
        self.assertTrue(witness["declared_theta_can_select_nonexact_sector"])
        self.assertFalse(witness["theta_selection_is_internal_source"])

    def test_source_matrix_blocks_false_passes(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        for candidate in (
            "positive_local_flatness_edge_action",
            "gauge_exact_boundary_phase_dynamics",
            "declared_theta_holonomy_source",
            "bridge_completed_sector_selector_theorem_target",
        ):
            self.assertIn(candidate, matrix)
            self.assertFalse(matrix[candidate]["passes_as_sector_selector_now"])
        self.assertIn("false pass", matrix["declared_theta_holonomy_source"]["verdict"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_sector_selector_candidates"], [])
        self.assertTrue(decision["nonexact_bit_carrier_still_available"])
        self.assertFalse(decision["local_even_class_selects_nonexact_sector"])
        self.assertTrue(decision["declared_theta_can_select_nonexact_sector"])
        self.assertFalse(decision["theta_selection_is_internal_source"])
        self.assertFalse(decision["boundary_phase_sector_selector_exported_now"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2664/S1614", MD.read_text(encoding="utf-8"))
        self.assertIn("P2664/S1614", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2664/S1614", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
