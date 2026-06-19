import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.py"
JSON_PATH = ROOT / "generated" / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.json"
MD_PATH = ROOT / "generated" / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.md"


class P2912GammaVariationalChainRuleSkeletonGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2911_input(self):
        self.assertEqual(self.payload["status"], "P2912_GAMMA_VARIATIONAL_CHAIN_RULE_SKELETON_GATE_READINESS_NO_EXPORT")
        self.assertTrue(self.acceptance["p2911_rechecked_localization_readiness"])

    def test_jacobian_counts(self):
        self.assertEqual(self.acceptance["edge_variable_count"], 144)
        self.assertEqual(self.acceptance["site_count"], 12)
        self.assertEqual(self.acceptance["jacobian_total_entries"], 1728)
        self.assertEqual(self.acceptance["jacobian_nonzero_entry_count"], 276)
        self.assertEqual(self.acceptance["jacobian_zero_entry_count"], 1452)
        self.assertEqual(len(self.objects["jacobian_nonzero_rows"]), 276)

    def test_target_edge_derivatives_and_covariance(self):
        self.assertEqual(self.acceptance["target_edge_0_5_nonzero_derivative_count"], 2)
        self.assertEqual(self.acceptance["target_edge_0_5_derivative_sites"], [0, 5])
        self.assertEqual(self.acceptance["target_edge_0_5_derivative_values"], ["1/2", "1/2"])
        self.assertEqual(self.acceptance["translation_covariance_failure_count"], 0)
        self.assertTrue(self.acceptance["finite_variational_chain_rule_skeleton_constructed"])

    def test_no_ltotal_or_eom_export(self):
        self.assertFalse(self.acceptance["strict_field_variable_provenance_exported"])
        self.assertFalse(self.acceptance["continuum_variational_chain_rule_exported"])
        self.assertFalse(self.acceptance["strict_gamma_9_5_source_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_variational_rule"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2912/S1862", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2912/S1862", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2912/S1862", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2912", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
