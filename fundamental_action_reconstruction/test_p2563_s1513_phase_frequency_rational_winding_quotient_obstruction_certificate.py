from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2563_s1513_phase_frequency_rational_winding_quotient_obstruction_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2563PhaseFrequencyRationalWindingQuotientObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_rational_winding_quotient_obstruction_certificate"]["theorem_export"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2563")
        self.assertEqual(self.payload["stage_id"], "S1513")
        self.assertIn("RATIONAL_WINDING_QUOTIENT_OBSTRUCTION", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2562_shortcut_nonpromotion_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])
        self.assertTrue(self.theorem["p2415_affine_nonautomorphism_inherited"])

    def test_exact_obstruction_and_bounded_search(self) -> None:
        self.assertEqual(len(self.theorem["exact_obstruction_rows"]), 2)
        self.assertTrue(self.theorem["exact_rational_winding_quotient_source_rejected"])
        for row in self.theorem["exact_obstruction_rows"]:
            self.assertFalse(row["exact_equality_possible_without_new_pi_cancelling_source"])
            self.assertIn("rational multiple of pi", row["proof_reason"])
            self.assertGreater(row["bounded_best_pi_multiple"]["abs_residual"], 0.0)
        search = self.theorem["bounded_pi_multiple_search"]
        self.assertEqual(search["numerator_bound"], 96)
        self.assertEqual(search["denominator_bound"], 96)

    def test_finite_phase_table_and_recommendation(self) -> None:
        self.assertEqual(self.theorem["finite_domain_size"], 12)
        self.assertEqual(len(self.theorem["finite_phase_table_d0_to_d11"]), 12)
        self.assertGreater(self.theorem["max_phase_gap_abs_d0_to_d11"], 0.0)
        self.assertGreater(self.theorem["max_cos_gap_abs_d0_to_d11"], 0.0)
        self.assertGreater(self.theorem["sign_mismatch_count_d0_to_d11"], 0)
        self.assertIn("non-winding strict phase/frequency source theorem", self.theorem["recommended_next_honest_step"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["rational_winding_quotient_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2563/S1513", MD.read_text(encoding="utf-8"))
        self.assertIn("P2563/S1513", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2563/S1513", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
