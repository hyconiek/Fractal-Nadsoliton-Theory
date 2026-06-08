from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2568_s1518_phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2568PhaseFrequencySemiboundedHessianRealizationNonuniquenessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate"]["theorem_export"]
        cls.audit = cls.theorem["semibounded_hessian_realization_audit"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2568")
        self.assertEqual(self.payload["stage_id"], "S1518")
        self.assertIn("SEMIBOUNDED_HESSIAN_REALIZATION_NONUNIQUENESS", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2567_minimal_saddle_obstruction_inherited"])
        self.assertTrue(self.theorem["p2566_stationarity_weight_cone_obstruction_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])

    def test_hessian_realization_nonuniqueness(self) -> None:
        self.assertEqual(self.audit["constraint_matrix_shape"], [5, 12])
        self.assertEqual(self.audit["constraint_matrix_rank"], 5)
        self.assertEqual(self.audit["solution_affine_nullity_for_fixed_hessian"], 7)
        self.assertTrue(self.audit["all_targets_realized_to_tolerance"])
        self.assertTrue(self.audit["both_local_max_and_local_min_hessians_realizable"])
        self.assertTrue(self.audit["all_realizations_use_signed_weights"])
        self.assertTrue(self.theorem["semibounded_hessian_can_be_realized_but_not_sourced"])
        self.assertTrue(self.theorem["hessian_sign_choice_is_extra_source_obligation"])
        for row in self.audit["target_rows"].values():
            self.assertLess(row["max_abs_residual"], 1e-10)
            self.assertTrue(row["has_positive_weights"])
            self.assertTrue(row["has_negative_weights"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["semibounded_hessian_weight_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2568/S1518", MD.read_text(encoding="utf-8"))
        self.assertIn("P2568/S1518", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2568/S1518", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
