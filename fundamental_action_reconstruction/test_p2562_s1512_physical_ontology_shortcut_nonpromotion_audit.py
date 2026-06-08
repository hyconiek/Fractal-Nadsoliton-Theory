from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2562_s1512_physical_ontology_shortcut_nonpromotion_audit import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2562PhysicalOntologyShortcutNonpromotionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["physical_ontology_shortcut_nonpromotion_audit"]["theorem_export"]
        cls.numeric = cls.theorem["numerical_shortcut_audit"]

    def test_identity_and_inherited_frontier(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2562")
        self.assertEqual(self.payload["stage_id"], "S1512")
        self.assertIn("SHORTCUT_NONPROMOTION", self.payload["status"])
        self.assertTrue(self.theorem["p2561_residual_two_key_frontier_inherited"])
        self.assertTrue(self.theorem["p2515_m2_signature_identified_not_sourced_inherited"])

    def test_heuristic_claims_not_promoted(self) -> None:
        self.assertEqual(self.theorem["heuristic_claim_count"], 3)
        self.assertEqual(self.theorem["promotable_heuristic_claim_count"], 0)
        self.assertTrue(self.theorem["physical_ontology_document_can_be_written_as_interpretation"])
        self.assertFalse(self.theorem["physical_ontology_document_can_replace_source_theorem"])
        self.assertFalse(self.theorem["heuristic_shortcut_changes_bridge_truth_table"])
        self.assertEqual(self.theorem["residual_atoms_before_shortcut"], self.theorem["residual_atoms_after_heuristic_shortcut"])

    def test_numerical_shortcut_audit(self) -> None:
        self.assertEqual(self.numeric["strict_parameters"], {"omega": "743/4000", "phi": "13/80", "eta": "9/5"})
        self.assertEqual(len(self.numeric["phase_rows_d1_to_d11"]), 11)
        self.assertAlmostEqual(self.numeric["eta_gap_strict_minus_legacy_linear_power"], 0.8, places=12)
        self.assertGreater(self.numeric["omega_gap_abs"], 0.0)
        self.assertGreater(self.numeric["phi_gap_abs"], 0.0)
        self.assertGreater(self.numeric["max_phase_gap_abs_d1_to_d11"], 0.0)
        self.assertGreater(self.numeric["max_cos_gap_abs_d1_to_d11"], 0.0)
        self.assertEqual(len(self.numeric["integer_rescaling_search_k1_to_k12"]), 12)
        self.assertTrue(self.numeric["no_exact_legacy_to_strict_numeric_identity"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["phenomenological_shortcut_promoted_to_source"])
        self.assertFalse(self.theorem["eta_fractal_tunneling_source_exported"])
        self.assertFalse(self.theorem["omega_nonlinear_resonance_source_exported"])
        self.assertFalse(self.theorem["hydrodynamic_m2_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2562/S1512", MD.read_text(encoding="utf-8"))
        self.assertIn("P2562/S1512", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2562/S1512", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
