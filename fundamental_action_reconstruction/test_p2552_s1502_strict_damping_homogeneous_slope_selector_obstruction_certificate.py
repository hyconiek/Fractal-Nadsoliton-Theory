from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2552_s1502_strict_damping_homogeneous_slope_selector_obstruction_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2552HomogeneousSlopeSelectorObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_homogeneous_slope_selector_obstruction_certificate"]["theorem_export"]
        cls.audit = cls.theorem["homogeneous_slope_selector_audit"]

    def test_identity_and_precursors(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2552")
        self.assertEqual(self.payload["stage_id"], "S1502")
        self.assertIn("HOMOGENEOUS_SLOPE_SELECTOR_OBSTRUCTION", self.payload["status"])
        self.assertTrue(self.theorem["p2543_slope_obstruction_inherited"])
        self.assertTrue(self.theorem["p2550_prime_log_basis_inherited"])
        self.assertTrue(self.theorem["p2551_post_prime_log_slope_obstruction_inherited"])

    def test_homogeneous_selector_dichotomy(self) -> None:
        self.assertEqual(self.audit["audited_row_count"], 6)
        self.assertEqual(self.audit["scale_invariant_rows"], ["ratio_edge_2_3", "ratio_edge_3_5", "ratio_edge_5_7", "ratio_edge_7_11"])
        self.assertEqual(self.audit["zero_selector_rows"], ["zero_prime_value_v2", "sum_prime_values_zero"])
        self.assertFalse(self.audit["any_row_uniquely_selects_strict_delta"])
        self.assertTrue(self.theorem["homogeneous_constraints_cannot_select_nonzero_strict_delta"])
        self.assertTrue(self.theorem["nonhomogeneous_anchor_required"])
        for row in self.audit["audited_rows"]:
            self.assertFalse(row["uniquely_selects_strict_delta_4_over_5"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertFalse(self.theorem["homogeneous_slope_selector_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("nonhomogeneous", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2552/S1502", MD.read_text(encoding="utf-8"))
        self.assertIn("P2552/S1502", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2552/S1502", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
