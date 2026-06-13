from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.py"
OUT = ROOT / "generated" / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.json"
MD = ROOT / "generated" / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.md"


class P2674TauSrcPair12TypedSeedO3AuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_o3_not_numbers_only(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "tau_src_anchor_content",
            "pair12_carrier_content",
            "chart_sensitive_atlas_content",
            "typed_seed_target_content",
            "actual_export_blocker_content",
            "collapse_blocker_content",
            "closure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2673_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2673_names_o3_missing"])
        self.assertTrue(upstream["p2673_distance_three"])
        self.assertTrue(upstream["p2673_no_boundary_phase_bit_target"])
        self.assertTrue(upstream["p2673_no_ltotal_reopening"])

    def test_o3_acceptance_matrix_separates_material_from_export(self) -> None:
        matrix = {item["id"]: item for item in self.payload["o3_acceptance_matrix"]}
        self.assertTrue(matrix["S1_tau_src_source_anchor_present"]["satisfied_now"])
        self.assertTrue(matrix["S2_surviving_pair12_carrier_present"]["satisfied_now"])
        self.assertTrue(matrix["S3_chart_sensitive_atlas_lane_present"]["satisfied_now"])
        self.assertFalse(matrix["S4_actual_chart_label_retaining_typed_seed_exported"]["satisfied_now"])
        self.assertFalse(matrix["S5_nonprojector_nonquotient_nonprelm_descent_law"]["satisfied_now"])
        self.assertFalse(matrix["S6_sigma_to_f301_typed_arrow_not_by_fiat"]["satisfied_now"])
        self.assertTrue(all(item["content_hits"] > 0 for item in matrix.values()))

    def test_finite_o3_lattice_blocks_false_pass(self) -> None:
        lattice = self.payload["finite_o3_lattice"]
        self.assertEqual(lattice["total_states_checked"], 64)
        self.assertEqual(lattice["passing_o3_states_count"], 1)
        self.assertFalse(lattice["current_state_passes_o3"])
        self.assertEqual(lattice["hamming_distance_from_o3_pass"], 3)
        self.assertEqual(
            set(lattice["missing_o3_subobligations_now"]),
            {
                "S4_actual_chart_label_retaining_typed_seed_exported",
                "S5_nonprojector_nonquotient_nonprelm_descent_law",
                "S6_sigma_to_f301_typed_arrow_not_by_fiat",
            },
        )

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertFalse(decision["current_state_passes_o3"])
        self.assertFalse(decision["o3_exported_now"])
        self.assertFalse(decision["boundary_square_arrow_allowed_next"])
        self.assertFalse(decision["sector_swap_invariant_allowed_next"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2674/S1624", MD.read_text(encoding="utf-8"))
        self.assertIn("P2674/S1624", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2674/S1624", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
