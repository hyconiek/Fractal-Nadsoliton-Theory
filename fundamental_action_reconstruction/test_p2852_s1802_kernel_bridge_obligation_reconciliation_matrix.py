import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.py"
JSON_PATH = ROOT / "generated" / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.json"
MD_PATH = ROOT / "generated" / "p2852_s1802_kernel_bridge_obligation_reconciliation_matrix.md"


class P2852KernelBridgeObligationReconciliationMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2852_KERNEL_BRIDGE_OBLIGATION_RECONCILIATION_MATRIX_NO_CLOSURE")
        inputs = self.payload["kernel_bridge_obligation_reconciliation"]["input_statuses_rechecked"]
        self.assertEqual(inputs["P2849"], "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")
        self.assertEqual(inputs["P2850"], "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE")
        self.assertEqual(inputs["P2851"], "P2851_AMPLITUDE_NORMALIZATION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")

    def test_only_syntax_target_satisfied(self):
        statuses = self.payload["kernel_bridge_obligation_reconciliation"]["target_statuses"]
        self.assertTrue(statuses["syntax_level_common_expression_basis"]["satisfied"])
        self.assertFalse(statuses["damping_compression_bridge_atom"]["satisfied"])
        self.assertFalse(statuses["amplitude_normalization_bridge_atom"]["satisfied"])
        self.assertFalse(statuses["full_kernel_completion_bridge"]["satisfied"])
        self.assertFalse(statuses["role_transfer"]["satisfied"])
        self.assertIn("eta_source_law_exported", statuses["damping_compression_bridge_atom"]["missing"])
        self.assertIn("phase_frequency_bridge_atom_exported", statuses["full_kernel_completion_bridge"]["missing"])

    def test_candidate_scoring_and_no_exports(self):
        scores = self.payload["kernel_bridge_obligation_reconciliation"]["candidate_next_atom_scores"]
        self.assertEqual(scores["phase_frequency_bridge_atom"]["replay_risk"], "low")
        self.assertGreater(scores["phase_frequency_bridge_atom"]["missing_obligation_reduction"], 0)
        self.assertEqual(scores["eml_syntax_replay"]["replay_risk"], "high")
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["role_transfer_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_bridge_reconciliation_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_full_kernel_bridge"])
        self.assertIn("P2852/S1802", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2852/S1802", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2852/S1802", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2852", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
