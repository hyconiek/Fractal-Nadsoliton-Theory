from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2596_s1546_nadsoliton_hydrodynamic_ir_m2_source_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2596NadsolitonHydrodynamicIRM2SourceTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_hydrodynamic_ir_m2_source_theorem"]["theorem_export"]
        cls.body = cls.theorem["nadsoliton_hydrodynamic_ir_m2_source_theorem"]

    def test_source_theorem_identity_and_export(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2596")
        self.assertEqual(self.payload["stage_id"], "S1546")
        self.assertIn("M2_SOURCE_THEOREM_EXPORTED", self.payload["status"])
        self.assertEqual(self.theorem["frontier_source_key_under_attack"], "m2_operator_signature_source")
        self.assertTrue(self.theorem["source_theorem_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])
        self.assertTrue(self.body["hydrodynamic_ir_rg_source_theorem_nonempty"])
        self.assertIn("operator of order m=2 is the unique IR selector", self.body["source_theorem_statement"])

    def test_hydrodynamic_rg_selection(self) -> None:
        self.assertEqual(self.body["hydrodynamic_fractal_dimension_Df"], "9/5")
        self.assertEqual(self.body["selected_operator_orders"], [2])
        self.assertTrue(self.body["unique_selected_order_is_m2"])
        rows = {row["operator_order_m"]: row for row in self.body["order_admissibility_rows"]}
        self.assertEqual(rows[0]["admissibility_status"], "forbidden")
        self.assertEqual(rows[1]["admissibility_status"], "forbidden")
        self.assertEqual(rows[2]["admissibility_status"], "selected")
        self.assertEqual(rows[4]["admissibility_status"], "irrelevant")
        self.assertTrue(all(value == "0" for value in self.body["higher_even_relative_ir_limits_vs_m2"].values()))

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2596/S1546", MD.read_text(encoding="utf-8"))
        self.assertIn("P2596/S1546", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2596/S1546", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
