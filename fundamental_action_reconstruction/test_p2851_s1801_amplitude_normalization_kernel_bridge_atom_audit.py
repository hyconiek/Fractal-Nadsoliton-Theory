import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.py"
JSON_PATH = ROOT / "generated" / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.json"
MD_PATH = ROOT / "generated" / "p2851_s1801_amplitude_normalization_kernel_bridge_atom_audit.md"


class P2851AmplitudeNormalizationKernelBridgeAtomAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2851_AMPLITUDE_NORMALIZATION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")
        audit = self.payload["amplitude_normalization_bridge_atom_audit"]
        self.assertEqual(audit["input_statuses_rechecked"]["P2850"], "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE")
        self.assertEqual(audit["parameters"]["alpha_geo_symbolic"], "4*ln(2)")

    def test_ratio_and_fit_reject_amplitude_only_bridge(self):
        audit = self.payload["amplitude_normalization_bridge_atom_audit"]
        self.assertFalse(audit["ratio_stats"]["constant_ratio"])
        self.assertTrue(audit["ratio_stats"]["sign_changes_present"])
        self.assertGreater(audit["least_squares_constant_amplitude_fit"]["max_abs_residual"], 1e-12)
        matrix = audit["two_point_constant_amplitude_matrix"]
        self.assertTrue(matrix["all_pairs_reject_constant_amplitude"])
        for row in matrix["rows"]:
            self.assertFalse(row["same_constant_amplitude"])

    def test_no_bridge_or_role_transfer_export(self):
        audit = self.payload["amplitude_normalization_bridge_atom_audit"]
        premise = audit["premise_matrix"]
        self.assertFalse(premise["accepted_as_amplitude_normalization_bridge_atom"])
        self.assertIn("alpha_geo_source_law_safe_for_strict_kernel", premise["missing_premises"])
        self.assertIn("role_transfer_theorem_available", premise["missing_premises"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["amplitude_normalization_bridge_atom_exported"])
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_amplitude_normalization_bridge_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_amplitude_normalization_bridge_atom"])
        self.assertIn("P2851/S1801", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2851/S1801", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2851/S1801", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2851", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
