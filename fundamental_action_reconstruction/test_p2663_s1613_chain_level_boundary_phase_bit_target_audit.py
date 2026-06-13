from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.py"
OUT = ROOT / "generated" / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json"
MD = ROOT / "generated" / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.md"


class P2663ChainLevelBoundaryPhaseBitTargetAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_and_upstream_consistency(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in ("chain_boundary_phase_content", "entropy_bit_target_content", "unit_map_content", "nonclosure_guard_content"):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2662_conditional_unique_scale_verified"])
        self.assertTrue(upstream["p2662_boundary_phase_bit_target_missing"])
        self.assertTrue(upstream["p2662_unconditional_uv_unit_not_selected"])
        self.assertTrue(upstream["p2662_qw2191_not_discharged"])

    def test_chain_level_witness_finds_carrier_but_not_unique_target(self) -> None:
        witness = self.payload["chain_level_boundary_phase_witness"]
        self.assertEqual(witness["all_z2_cochain_count"], 64)
        self.assertEqual(witness["triangle_flat_cochain_count"], 32)
        self.assertEqual(witness["exact_coboundary_square_holonomy_values"], [0])
        self.assertEqual(witness["square_holonomy_values_after_flatness"], [0, 1])
        self.assertTrue(witness["nonexact_sector_needed_for_nonzero_bit"])
        self.assertFalse(witness["chain_law_derives_unique_n_bits_without_sector_choice"])
        self.assertTrue(witness["boundary_of_boundary_zero_gives_only_flatness_not_entropy_level"])

    def test_source_matrix_blocks_false_passes(self) -> None:
        matrix = {row["candidate"]: row for row in self.payload["source_candidate_matrix"]}
        for candidate in ("boundary_of_boundary_zero_as_bit_target", "exact_coboundary_phase_law", "nonexact_square_holonomy_bit", "declared_cycle_basis_entropy_target"):
            self.assertIn(candidate, matrix)
            self.assertFalse(matrix[candidate]["passes_as_boundary_phase_bit_target_now"])
        self.assertIn("promising carrier", matrix["nonexact_square_holonomy_bit"]["verdict"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_boundary_phase_bit_target_candidates"], [])
        self.assertTrue(decision["nonexact_bit_carrier_found"])
        self.assertFalse(decision["unique_n_bits_derived_without_sector_choice"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2663/S1613", MD.read_text(encoding="utf-8"))
        self.assertIn("P2663/S1613", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2663/S1613", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
