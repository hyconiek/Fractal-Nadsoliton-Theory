import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2853_s1803_phase_frequency_bridge_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2853_s1803_phase_frequency_bridge_source_audit.json"
MD_PATH = ROOT / "generated" / "p2853_s1803_phase_frequency_bridge_source_audit.md"


class P2853PhaseFrequencyBridgeSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["phase_frequency_bridge_source_audit"]
        cls.summary = cls.audit["summary"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2853_PHASE_FREQUENCY_BRIDGE_SOURCE_AUDIT_NO_CLOSURE")
        inputs = self.audit["input_statuses_rechecked"]
        self.assertEqual(inputs["P2852"], "P2852_KERNEL_BRIDGE_OBLIGATION_RECONCILIATION_MATRIX_NO_CLOSURE")

    def test_finite_transport_positive_but_not_source(self):
        self.assertTrue(self.summary["continuous_affine_phase_transport_exact"])
        self.assertLessEqual(self.summary["max_abs_affine_transport_residual"], 1e-12)
        self.assertTrue(self.summary["phase_factor_bits_nonconstant"])
        self.assertTrue(self.payload["acceptance_matrix"]["exports_phase_frequency_transport_witness"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_phase_frequency_bridge_source_atom"])

    def test_rejects_discrete_and_scalar_shortcuts(self):
        self.assertTrue(self.summary["no_z12_unit_offset_reindex_matches_strict_sign_pattern"])
        self.assertEqual(self.summary["best_z12_unit_offset_mismatch_count"], 2)
        self.assertTrue(self.summary["scalar_phase_replacement_fails"])
        self.assertGreater(self.summary["scalar_phase_best_fit"]["max_abs_residual"], 0.1)
        self.assertTrue(self.summary["same_integer_coordinate_phase_identity_fails"])
        self.assertTrue(self.summary["constant_phase_shift_only_fails"])

    def test_candidate_matrix_has_no_export(self):
        self.assertEqual(self.audit["accepted_candidate_count"], 0)
        candidates = {row["candidate"]: row for row in self.audit["source_candidate_matrix"]}
        self.assertTrue(candidates["continuous_affine_phase_coordinate_transport"]["finite_witness_passes"])
        self.assertFalse(candidates["continuous_affine_phase_coordinate_transport"]["exports_strict_source_law"])
        self.assertFalse(candidates["z12_unit_offset_reindexing"]["finite_witness_passes"])
        self.assertFalse(candidates["phase_sign_gf2_bookkeeping"]["exports_strict_source_law"])

    def test_no_false_closure_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_phase_frequency_bridge_source_obstruction_audit"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["phase_frequency_bridge_source_atom_exported"])
        self.assertFalse(flags["strict_omega_phi_source_law_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["role_transfer_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])
        self.assertIn("P2853/S1803", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2853/S1803", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2853/S1803", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2853", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
