from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2555_s1505_legacy_to_strict_damping_denominator_nonrenormalization_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2555LegacyToStrictDampingDenominatorNonrenormalizationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_denominator_nonrenormalization_certificate"]["theorem_export"]

    def test_identity_and_precursor(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2555")
        self.assertEqual(self.payload["stage_id"], "S1505")
        self.assertIn("DAMPING_DENOMINATOR_NONRENORMALIZATION", self.payload["status"])
        self.assertTrue(self.theorem["p2554_bridge_reorientation_inherited"])

    def test_linear_vs_nonlinear_denominator_audit(self) -> None:
        self.assertEqual(self.theorem["beta_tors_candidates_audited"], [0.01, 0.05])
        self.assertTrue(self.theorem["legacy_linear_second_difference_zero_for_all_candidates"])
        self.assertTrue(self.theorem["strict_nonlinear_second_difference_positive"])
        self.assertTrue(self.theorem["raw_denominator_identity_refuted_on_domain"])
        self.assertTrue(self.theorem["constant_amplitude_absorption_refuted_on_domain"])
        self.assertTrue(self.theorem["beta_tors_to_beta_eta_scalar_renormalization_refuted_for_audited_class"])
        for row in self.theorem["damping_denominator_nonrenormalization_rows"]:
            self.assertLess(row["legacy_second_difference_max_abs"], 1e-12)
            self.assertGreater(row["strict_second_difference_min"], 0.0)
            self.assertFalse(row["constant_amplitude_can_absorb_denominator_mismatch"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertFalse(self.theorem["beta_tors_to_beta_eta_translation_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["legacy_to_strict_completion_bridge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_claimed"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("non-scalar nonlinear source", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2555/S1505", MD.read_text(encoding="utf-8"))
        self.assertIn("P2555/S1505", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2555/S1505", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
