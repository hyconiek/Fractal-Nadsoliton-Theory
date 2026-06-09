from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2615_s1565_p2605_linear_slice_nonbridge_obstruction import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2615P2605LinearSliceNonbridgeObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2605_linear_slice_nonbridge_obstruction"]["theorem_export"]

    def test_identity_and_proof_shape(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2615")
        self.assertEqual(self.payload["stage_id"], "S1565")
        self.assertIn("NONBRIDGE_OBSTRUCTION", self.payload["status"])
        proof_text = " ".join(self.theorem["exact_two_node_obstruction"]["proof_steps"])
        self.assertIn("eta=1", proof_text)
        self.assertIn("beta=beta_tors", proof_text)

    def test_nonlinear_candidates_fail_two_node_constant_beta_completion(self) -> None:
        rows = self.theorem["nonlinear_candidate_rows"]
        self.assertGreaterEqual(len(rows), 2)
        for row in rows:
            self.assertFalse(row["eta_equals_one"])
            self.assertFalse(row["accepts_single_constant_beta_completion_on_two_nodes"])
            self.assertGreater(row["two_anchor_beta_mismatch"], 0.0)
            self.assertGreater(row["max_abs_residual_when_beta_fit_at_d1"], 0.0)

    def test_quarantine_and_retention_logic(self) -> None:
        self.assertTrue(self.theorem["p2605_quarantine_before_p2615"])
        self.assertTrue(self.theorem["p2605_linear_slice_reclassified_as_boundary_negative_control"])
        self.assertTrue(self.theorem["p2605_quarantine_retained_by_p2615"])
        self.assertTrue(self.theorem["p2606_nonlinear_residual_component_retained"])
        self.assertTrue(self.theorem["p2601_p2602_revalidation_inherited"])
        self.assertIn("P2605", self.theorem["remaining_p2610_quarantines_after_p2615"])

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["p2605_full_bridge_revalidated"])
        self.assertFalse(self.theorem["p2607_bridge_completion_revalidated"])
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2615/S1565", MD.read_text(encoding="utf-8"))
        self.assertIn("P2615/S1565", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2615/S1565", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
