from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.py"
OUT = ROOT / "generated" / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json"
MD = ROOT / "generated" / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.md"


class P2689EntropyReferenceCellBitToLengthUvUnitObligationMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_grep_covers_entropy_uv_unit_frontier(self) -> None:
        audit = self.payload["content_grep"]
        self.assertIn("entropy reference-cell/bit-to-length", audit["mode"])
        for key in (
            "p2688_selected_p2689",
            "p2662_conditional_unit_map",
            "p2663_one_bit_carrier",
            "scale_orbit_beta_source",
            "forbidden_replays",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_state_reads_use_p2688_p2662_p2663_and_scale_orbit(self) -> None:
        state = self.payload["state_reads"]
        self.assertTrue(state["p2688_selected_p2689"])
        self.assertTrue(state["p2662_conditional_unique_scale_verified"])
        self.assertFalse(state["p2662_unconditional_uv_unit_selected"])
        self.assertIn("intrinsic_pre_normalization_entropy_measure", state["p2662_missing_premises"])
        self.assertTrue(state["p2663_nonexact_sector_needed"])
        self.assertFalse(state["p2663_unique_bits_without_sector_choice"])
        self.assertEqual(state["p2663_exact_coboundary_values"], [0])
        self.assertFalse(state["p2650_any_candidate_passes"])
        self.assertFalse(state["p2653_canonical_unit_exported"])
        self.assertTrue(state["p2653_scale_orbit_equivalence_verified"])

    def test_symbolic_unit_map_is_conditional_positive_scale(self) -> None:
        sym = self.payload["symbolic_unit_map_check"]
        self.assertEqual(sym["entropy_law"], "Df*log(a) + H0")
        self.assertEqual(sym["bit_target"], "N*log(2)")
        self.assertIn("exp", sym["selected_scale_formula"])
        self.assertTrue(sym["unique_positive_scale_if_premises_supplied"])
        self.assertTrue(sym["conditional_math_passes"])
        self.assertTrue(all(row["positive"] for row in sym["sample_scales"]))

    def test_obligation_matrix_blocks_unconditional_uv_unit(self) -> None:
        obligations = self.payload["obligation_matrix"]
        self.assertEqual(len(obligations), 5)
        self.assertTrue(all(not row["satisfied_now"] for row in obligations))
        ids = {row["obligation_id"] for row in obligations}
        self.assertIn("intrinsic_entropy_reference_cell_or_entropy_zero", ids)
        self.assertIn("selector_free_boundary_phase_bit_target", ids)
        self.assertIn("bit_to_length_or_action_unit_map", ids)
        self.assertIn("scale_orbit_quotient_single_uv_unit", ids)
        self.assertIn("target_independent_beta_or_z_beta_source", ids)

    def test_decision_is_bounded_no_go_without_exports_and_updates_docs(self) -> None:
        decision = self.payload["decision"]
        self.assertFalse(decision["unconditional_uv_unit_exported"])
        self.assertTrue(decision["bounded_no_go_for_current_canonical_unit_atom"])
        self.assertFalse(decision["beta_or_z_beta_source_exported_now"])
        self.assertFalse(decision["role_transfer_started_now"])
        self.assertFalse(decision["toe_closed_now"])
        self.assertIn("P2690", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2689/S1639", MD.read_text(encoding="utf-8"))
        self.assertIn("P2689/S1639", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2689/S1639", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2689/S1639", (ROOT.parent / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
