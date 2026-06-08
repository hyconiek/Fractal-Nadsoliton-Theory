from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2564_s1514_phase_frequency_finite_sign_cell_nonidentifiability_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2564PhaseFrequencyFiniteSignCellNonidentifiabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["phase_frequency_finite_sign_cell_nonidentifiability_certificate"]["theorem_export"]

    def test_identity_and_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2564")
        self.assertEqual(self.payload["stage_id"], "S1514")
        self.assertIn("FINITE_SIGN_CELL_NONIDENTIFIABILITY", self.payload["status"])
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_phase_frequency_source")
        self.assertTrue(self.theorem["p2563_rational_winding_obstruction_inherited"])
        self.assertTrue(self.theorem["p2561_phase_frequency_residual_atom_inherited"])
        self.assertTrue(self.theorem["p2415_affine_nonautomorphism_inherited"])

    def test_open_sign_cell_and_grid_witnesses(self) -> None:
        self.assertEqual(len(self.theorem["strict_phase_rows_d0_to_d11"]), 12)
        self.assertEqual(len(self.theorem["strict_sign_pattern_d0_to_d11"]), 12)
        box = self.theorem["certified_open_sign_cell_box"]
        self.assertGreater(box["epsilon_omega"], 0.0)
        self.assertGreater(box["epsilon_phi"], 0.0)
        self.assertGreater(box["box_area"], 0.0)
        self.assertEqual(self.theorem["grid_witness_count"], 25)
        self.assertEqual(self.theorem["grid_witnesses_preserving_sign_count"], 25)
        self.assertGreater(self.theorem["nontrivial_same_sign_witness_count"], 0)
        self.assertFalse(self.theorem["finite_sign_pattern_selects_unique_phase_frequency_tuple"])
        self.assertTrue(self.theorem["finite_sign_pattern_has_open_continuum_of_phase_frequency_realizations"])

    def test_negative_controls_and_recommendation(self) -> None:
        self.assertIn("stronger selector candidate", self.theorem["recommended_next_honest_step"])
        self.assertFalse(self.theorem["finite_sign_pattern_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_docs(self) -> None:
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2564/S1514", MD.read_text(encoding="utf-8"))
        self.assertIn("P2564/S1514", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2564/S1514", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
