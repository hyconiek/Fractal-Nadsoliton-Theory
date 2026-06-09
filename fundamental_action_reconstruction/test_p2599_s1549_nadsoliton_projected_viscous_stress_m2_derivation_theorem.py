from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2599_s1549_nadsoliton_projected_viscous_stress_m2_derivation_theorem import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2599NadsolitonProjectedViscousStressM2DerivationTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["nadsoliton_projected_viscous_stress_m2_derivation_theorem"]["theorem_export"]
        cls.body = cls.theorem["nadsoliton_projected_viscous_stress_m2_derivation_theorem"]

    def test_identity_and_source_export(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2599")
        self.assertEqual(self.payload["stage_id"], "S1549")
        self.assertIn("PROJECTED_VISCOUS_STRESS_M2_DERIVATION", self.payload["status"])
        self.assertTrue(self.theorem["p2598_locality_source_theorem_inherited"])
        self.assertTrue(self.theorem["source_theorem_exported"])
        self.assertTrue(self.theorem["m2_operator_signature_source_exported"])
        self.assertIn("Fourier generator on transverse modes is exactly -mu |k|^2", self.body["source_theorem_statement"])

    def test_projected_viscous_symbol(self) -> None:
        self.assertEqual(self.body["selected_operator_order"], 2)
        self.assertEqual(self.body["symbolic_stress_derivation"]["operator_order_from_symbol"], 2)
        self.assertTrue(self.body["all_projectors_idempotent"])
        self.assertTrue(self.body["all_longitudinal_modes_removed"])
        self.assertTrue(self.body["all_transverse_modes_have_minus_mu_k2_symbol"])
        self.assertEqual(len(self.body["projected_viscous_sample_rows"]), 4)
        for row in self.body["projected_viscous_sample_rows"]:
            self.assertEqual(row["transverse_basis_count"], 2)
            self.assertTrue(row["all_transverse_checks_zero"])

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["fractional_laplacian_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_theorem"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2599/S1549", MD.read_text(encoding="utf-8"))
        self.assertIn("P2599/S1549", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2599/S1549", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
