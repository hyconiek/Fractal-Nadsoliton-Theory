from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2617_s1567_p2606_exponent_semantics_reclassification import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2617P2606ExponentSemanticsReclassificationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2606_exponent_semantics_reclassification"]["theorem_export"]

    def test_identity_and_semantics(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2617")
        self.assertEqual(self.payload["stage_id"], "S1567")
        self.assertIn("EXPONENT_SEMANTICS_RECLASSIFIED", self.payload["status"])
        self.assertEqual(self.theorem["delta_codimension"], "4/5")
        self.assertEqual(self.theorem["eta_strict_denominator"], "9/5")
        self.assertTrue(self.theorem["eta_strict_equals_one_plus_delta"])

    def test_p2606_quarantined_as_strict_residual_but_probe_retained(self) -> None:
        self.assertTrue(self.theorem["p2606_used_delta_as_denominator_eta"])
        self.assertGreater(self.theorem["max_abs_strict_eta_9_5_minus_p2606_eta_4_5_kernel"], 1e-6)
        self.assertTrue(self.theorem["p2606_strict_residual_interpretation_quarantined_by_p2617"])
        self.assertTrue(self.theorem["p2606_codimension_probe_retained"])

    def test_proof_mentions_correct_relation(self) -> None:
        proof_text = " ".join(self.theorem["exact_semantics_proof"]["proof_steps"])
        self.assertIn("delta = D_f - 1 = 4/5", proof_text)
        self.assertIn("eta_strict = 1 + delta", proof_text)
        self.assertIn("9/5", proof_text)

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["p2606_strict_side_residual_addition_revalidated"])
        self.assertFalse(self.theorem["p2607_bridge_completion_revalidated"])
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2617/S1567", MD.read_text(encoding="utf-8"))
        self.assertIn("P2617/S1567", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2617/S1567", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
