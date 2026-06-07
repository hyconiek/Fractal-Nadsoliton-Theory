from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2561PostDampingResidualBridgeTwoKeyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_completion_post_damping_residual_bridge_two_key_certificate"]["theorem_export"]

    def test_residual_truth_table(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2561")
        self.assertEqual(self.payload["stage_id"], "S1511")
        self.assertTrue(self.theorem["p2502_unique_bridge_triple_inherited"])
        self.assertTrue(self.theorem["p2560_constant_log_source_obstruction_inherited"])
        self.assertEqual(self.theorem["current_bridge_accepting_rows"], 0)
        self.assertEqual(self.theorem["best_case_damping_truth_table_row_count"], 4)
        self.assertEqual(self.theorem["best_case_damping_bridge_accepting_row_count"], 1)
        self.assertEqual(self.theorem["residual_atoms_after_hypothetical_damping_source"], ["strict_dynamical_source_for_A_P_D", "strict_phase_frequency_source"])
        self.assertTrue(self.theorem["residual_atoms_are_jointly_required_after_damping"])
        self.assertEqual(len(self.theorem["single_residual_omission_witnesses"]), 2)

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("phase/frequency/topological-bit passage", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2561/S1511", MD.read_text(encoding="utf-8"))
        self.assertIn("P2561/S1511", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2561/S1511", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
