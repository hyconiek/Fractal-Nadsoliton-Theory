from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2597_s1547_nadsoliton_hydrodynamic_rg_m2_robustness_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2597NadsolitonHydrodynamicRGM2RobustnessTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_hydrodynamic_rg_m2_robustness_theorem"]["theorem_export"]
        cls.body = cls.theorem["nadsoliton_hydrodynamic_rg_m2_robustness_theorem"]

    def test_identity_and_source_export(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2597")
        self.assertEqual(self.payload["stage_id"], "S1547")
        self.assertIn("M2_ROBUSTNESS_THEOREM_EXPORTED", self.payload["status"])
        self.assertTrue(self.theorem["p2596_m2_source_theorem_inherited"])
        self.assertTrue(self.theorem["source_theorem_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])
        self.assertIn("RG-robust", self.body["robust_source_theorem_statement"])

    def test_rg_robustness_grid(self) -> None:
        self.assertEqual(self.body["D_f_interval_audited"], ["17/10", "19/10"])
        self.assertEqual(self.body["selected_operator_orders"], [2])
        self.assertTrue(self.body["unique_selected_order_is_m2_for_all_Df_grid"])
        self.assertTrue(self.body["higher_even_orders_all_have_negative_relative_eigenvalue"])
        self.assertTrue(self.body["m0_m1_and_odd_orders_forbidden_for_all_Df_grid"])
        self.assertTrue(self.body["relative_rg_eigenvalue_independent_of_Df"])
        self.assertEqual(self.body["minimum_ir_relevance_gap_to_next_even_operator"], 2)
        rows_m4 = [row for row in self.body["rg_rows"] if row["operator_order_m"] == 4]
        self.assertTrue(all(row["relative_rg_eigenvalue_y_m_minus_y_2"] == "-2" for row in rows_m4))

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2597/S1547", MD.read_text(encoding="utf-8"))
        self.assertIn("P2597/S1547", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2597/S1547", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
