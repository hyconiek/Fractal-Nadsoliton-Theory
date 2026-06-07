from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2553_s1503_strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2553NonhomogeneousAnchorConstantEquivalenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["strict_damping_nonhomogeneous_anchor_constant_equivalence_certificate"]["theorem_export"]
        cls.audit = cls.theorem["nonhomogeneous_anchor_constant_equivalence_audit"]

    def test_identity_and_precursors(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2553")
        self.assertEqual(self.payload["stage_id"], "S1503")
        self.assertIn("NONHOMOGENEOUS_ANCHOR_CONSTANT_EQUIVALENCE", self.payload["status"])
        self.assertTrue(self.theorem["p2528_prime_anchor_equivalence_inherited"])
        self.assertTrue(self.theorem["p2551_post_prime_log_slope_obstruction_inherited"])
        self.assertTrue(self.theorem["p2552_homogeneous_selector_obstruction_inherited"])

    def test_anchor_constant_equivalence(self) -> None:
        self.assertEqual(self.audit["anchor_row_count"], 7)
        self.assertTrue(self.audit["all_anchor_rows_have_nonzero_log_dot"])
        for row in self.audit["anchor_rows"]:
            self.assertEqual(row["selector_formula_if_k_supplied"], "delta = k/(c·log(p))")
            self.assertEqual(row["delta_selected_by_strict_constant"], "4/5")
            self.assertTrue(row["constant_is_equivalent_to_slope_value"])
        self.assertTrue(self.theorem["nonhomogeneous_anchor_reduces_to_constant_source_obligation"])

    def test_misanchor_witnesses_and_negative_controls(self) -> None:
        self.assertEqual(len(self.audit["constant_misanchor_witnesses"]), 21)
        self.assertTrue(all(not witness["selects_strict_delta_4_over_5"] for witness in self.audit["constant_misanchor_witnesses"] if witness["constant_label"] != "strict_constant"))
        self.assertFalse(self.theorem["nonhomogeneous_anchor_constant_source_exported"])
        self.assertFalse(self.theorem["slope_value_or_prime_anchor_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("independently sourced", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2553/S1503", MD.read_text(encoding="utf-8"))
        self.assertIn("P2553/S1503", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2553/S1503", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
