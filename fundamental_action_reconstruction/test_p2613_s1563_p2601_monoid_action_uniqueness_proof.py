from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2613_s1563_p2601_monoid_action_uniqueness_proof import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2613P2601MonoidActionUniquenessProofTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2601_monoid_action_uniqueness_proof"]["theorem_export"]
        cls.audit = cls.theorem["finite_consistency_audit"]

    def test_identity_and_prompt_verdict(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2613")
        self.assertEqual(self.payload["stage_id"], "S1563")
        self.assertIn("STRICT_MONOID_ACTION_UNIQUENESS", self.payload["status"])
        self.assertIn("ACCEPTED_WITH_SCOPE", self.theorem["prompt_professorial_verdict"])
        self.assertEqual(len(self.theorem["monoid_action_axioms"]), 5)

    def test_y1_zero_proved_and_p2601_lifted(self) -> None:
        self.assertTrue(self.theorem["p2601_quarantine_before_p2613"])
        self.assertTrue(self.theorem["p2601_monoid_action_uniqueness_revalidated"])
        self.assertTrue(self.theorem["unital_normalization_y1_zero_source_exported"])
        self.assertTrue(self.theorem["multiplicative_character_law_source_revalidated"])
        self.assertTrue(self.theorem["p2601_quarantine_lifted_by_p2613"])
        self.assertNotIn("P2601", self.theorem["remaining_p2610_quarantines_after_p2613"])

    def test_algebraic_and_computational_checks(self) -> None:
        self.assertIn("y_1 + y_1", " ".join(self.theorem["proof_core"]["algebraic_steps"]))
        self.assertIn("D_f does not enter", self.theorem["proof_core"]["df_independence"])
        self.assertTrue(self.audit["only_zero_candidate_passes_additive_identity"])
        self.assertTrue(self.audit["only_zero_candidate_passes_attenuation_identity"])
        self.assertGreater(self.audit["product_rows_count"], 12)

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["p2602_prime_log_source_revalidated"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_fully_revalidated"])
        self.assertFalse(self.theorem["bridge_completion_revalidated"])
        self.assertFalse(self.theorem["role_bearing_ltotal_reenabled"])
        self.assertTrue(self.theorem["p2612_gf2_obstruction_respected"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2613/S1563", MD.read_text(encoding="utf-8"))
        self.assertIn("P2613/S1563", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2613/S1563", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
