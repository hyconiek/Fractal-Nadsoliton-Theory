from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.py"
OUT = ROOT / "generated" / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.json"
MD = ROOT / "generated" / "p2667_s1617_pair12_boundary_orientation_reversal_no_go_audit.md"


class P2667Pair12BoundaryOrientationReversalNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_orientation_and_boundary_cycle(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "orientation_convention_content",
            "pair12_positive_branch_content",
            "boundary_cycle_functor_content",
            "symmetry_reversal_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2666_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2666_pair12_split_exported"])
        self.assertTrue(upstream["p2666_pair2_positive"])
        self.assertTrue(upstream["p2666_possible_sector_one_mapping_exists"])
        self.assertTrue(upstream["p2666_no_passing_descent_mapping"])
        self.assertTrue(upstream["p2666_strict_current_descent_not_exported"])

    def test_orientation_reversal_witness_blocks_conventional_map(self) -> None:
        witness = self.payload["orientation_reversal_witness"]
        self.assertEqual(witness["candidate_mapping_count"], 2)
        self.assertTrue(witness["sector_one_mapping_exists"])
        self.assertTrue(witness["all_reversal_pairs_unforbidden"])
        self.assertFalse(witness["canonical_orientation_map_exported_now"])
        self.assertEqual(len(witness["candidate_mappings"]), 2)
        self.assertTrue(any(row["selects_nonexact_sector"] for row in witness["candidate_mappings"]))
        self.assertTrue(all(not row["cross_invariant_pair2_equals_sector1_exported"] for row in witness["candidate_mappings"]))

    def test_source_matrix_and_nonclosure(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        for candidate in (
            "pair2_positive_label_as_boundary_sector_1",
            "boundary_sector_label_1_as_intrinsic_orientation",
            "sector_swap_reversal_symmetry",
            "future_cross_invariant_pair2_to_sector1",
        ):
            self.assertIn(candidate, matrix)
            self.assertFalse(matrix[candidate]["passes_now"])
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_orientation_source_candidates"], [])
        self.assertFalse(decision["canonical_orientation_map_exported_now"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2667/S1617", MD.read_text(encoding="utf-8"))
        self.assertIn("P2667/S1617", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2667/S1617", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
