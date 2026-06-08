from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2598_s1548_nadsoliton_hydrodynamic_locality_fractional_competitor_exclusion_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2598NadsolitonHydrodynamicLocalityFractionalCompetitorExclusionTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_hydrodynamic_locality_fractional_competitor_exclusion_theorem"]["theorem_export"]
        cls.body = cls.theorem["nadsoliton_hydrodynamic_locality_fractional_competitor_exclusion_theorem"]

    def test_identity_and_source_export(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2598")
        self.assertEqual(self.payload["stage_id"], "S1548")
        self.assertIn("FRACTIONAL_COMPETITOR_EXCLUSION_THEOREM_EXPORTED", self.payload["status"])
        self.assertTrue(self.theorem["p2597_rg_m2_robustness_inherited"])
        self.assertTrue(self.theorem["source_theorem_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])
        self.assertIn("fractional competitors alpha<2 are excluded", self.body["source_theorem_statement"])

    def test_fractional_competitors_excluded(self) -> None:
        self.assertEqual(self.body["selected_orders"], ["2"])
        self.assertTrue(self.body["unique_selected_order_is_2"])
        self.assertTrue(self.body["all_fractional_alpha_below_2_excluded_by_locality"])
        self.assertTrue(self.body["all_alpha_above_2_irrelevant"])
        rows = {row["operator_order_alpha"]: row for row in self.body["fractional_order_rows"]}
        self.assertEqual(rows["1/2"]["admissibility_status"], "forbidden_nonlocal_fractional_competitor")
        self.assertEqual(rows["3/2"]["admissibility_status"], "forbidden_nonlocal_fractional_competitor")
        self.assertEqual(rows["2"]["admissibility_status"], "selected_local_laplacian")
        self.assertEqual(rows["5/2"]["admissibility_status"], "irrelevant_higher_local_or_hyperlocal_order")

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["fractional_laplacian_source_exported"])
        self.assertFalse(self.theorem["nonlocal_stable_transport_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2598/S1548", MD.read_text(encoding="utf-8"))
        self.assertIn("P2598/S1548", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2598/S1548", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
