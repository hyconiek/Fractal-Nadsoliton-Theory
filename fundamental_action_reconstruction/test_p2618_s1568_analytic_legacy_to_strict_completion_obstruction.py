from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2618_s1568_analytic_legacy_to_strict_completion_obstruction import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2618AnalyticLegacyToStrictCompletionObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["analytic_legacy_to_strict_completion_obstruction"]["theorem_export"]

    def test_identity_and_verdict(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2618")
        self.assertEqual(self.payload["stage_id"], "S1568")
        self.assertIn("ANALYTIC_COMPLETION_MAP_PARTIAL", self.payload["status"])
        self.assertEqual(self.theorem["result_type"], "fundamental_obstruction_to_full_a_priori_completion_map_under_current_sources")
        self.assertFalse(self.theorem["full_completion_map_exported"])

    def test_damping_exponent_partial_result_and_scalar_obstruction(self) -> None:
        damping = self.theorem["damping_partial_result"]
        self.assertEqual(damping["df_used"], "9/5")
        self.assertEqual(damping["delta_codimension"], "4/5")
        self.assertEqual(damping["eta_strict_denominator"], "9/5")
        self.assertEqual(damping["eta_relation"], "eta_strict = 1 + delta = D_f")
        proof_text = " ".join(self.theorem["damping_nonrenormalization_proof"]["proof_steps"])
        self.assertIn("c=1", proof_text)
        self.assertIn("d^(4/5)", proof_text)
        self.assertIn("impossible", proof_text)
        self.assertTrue(all(row["strict_exceeds_legacy_for_beta_1_beta_tors_0_01"] for row in self.theorem["damping_sample_sanity_check"]))

    def test_phase_selector_obstruction(self) -> None:
        selector = self.theorem["phase_selector_obstruction_proof"]
        proof_text = " ".join(selector["proof_steps"])
        self.assertIn("sigma = -sigma", proof_text)
        self.assertIn("impossible", proof_text)
        for row in selector["orientation_reversal_sanity_check"]:
            self.assertFalse(row["orientation_invariant_selector_condition_s_equals_reversed_s"])
            self.assertFalse(row["odd_phase_sign_survives_without_orientation_source"])
        self.assertFalse(self.theorem["gf2_free_strict_selector_exported"])

    def test_scope_guards(self) -> None:
        self.assertFalse(self.theorem["beta_tors_to_beta_theorem_exported"])
        self.assertFalse(self.theorem["p2607_bridge_completion_revalidated"])
        self.assertFalse(self.theorem["p2608_role_bearing_ltotal_reenabled"])
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_packet"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertFalse(self.theorem["apd_source_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2618/S1568", MD.read_text(encoding="utf-8"))
        self.assertIn("P2618/S1568", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2618/S1568", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
