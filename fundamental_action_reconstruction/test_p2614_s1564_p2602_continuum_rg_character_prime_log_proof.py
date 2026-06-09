from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2614_s1564_p2602_continuum_rg_character_prime_log_proof import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2614P2602ContinuumRgCharacterPrimeLogProofTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["p2602_continuum_rg_character_prime_log_proof"]["theorem_export"]

    def test_identity_and_axioms(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2614")
        self.assertEqual(self.payload["stage_id"], "S1564")
        self.assertIn("CONTINUUM_RG_CHARACTER", self.payload["status"])
        self.assertEqual(len(self.theorem["continuum_rg_character_axioms"]), 3)
        self.assertIn("Cauchy", " ".join(self.theorem["proof_core"]["algebraic_steps"]))

    def test_prime_log_revalidated(self) -> None:
        self.assertTrue(self.theorem["p2602_quarantine_before_p2614"])
        self.assertTrue(self.theorem["p2601_revalidated_inherited"])
        self.assertTrue(self.theorem["p2602_prime_log_proportionality_revalidated"])
        self.assertTrue(self.theorem["p2602_quarantine_lifted_by_p2614"])
        self.assertNotIn("P2602", self.theorem["remaining_p2610_quarantines_after_p2614"])

    def test_continuum_vs_discrete_counterexamples(self) -> None:
        self.assertTrue(all(row["prime_log_proportionality_accepts"] for row in self.theorem["log_character_gamma_rows"]))
        for row in self.theorem["noncontinuum_prime_character_counterexamples"]:
            self.assertTrue(row["passes_integer_multiplicative_character"])
            self.assertFalse(row["passes_continuum_log_character"])
            self.assertGreater(row["max_continuum_log_embedding_residual"], 0.0)

    def test_scope_guards(self) -> None:
        self.assertTrue(self.theorem["p2603_conditional_slope_retained"])
        self.assertTrue(self.theorem["strict_damping_beta_eta_source_revalidated_under_df_9_5_scope"])
        self.assertFalse(self.theorem["p2607_bridge_completion_revalidated"])
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2614/S1564", MD.read_text(encoding="utf-8"))
        self.assertIn("P2614/S1564", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2614/S1564", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
