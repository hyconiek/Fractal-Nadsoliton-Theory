import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.py"
JSON_PATH = ROOT / "generated" / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.json"
MD_PATH = ROOT / "generated" / "p2849_s1799_damping_compression_kernel_bridge_atom_audit.md"


class P2849DampingCompressionKernelBridgeAtomAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")
        audit = self.payload["damping_compression_bridge_atom_audit"]
        self.assertEqual(audit["input_statuses_rechecked"]["P2848"], "P2848_COUPLING_COEFFICIENT_UNIT_SOURCE_LAW_AUDIT_NO_CLOSURE")
        self.assertEqual(audit["parameters"]["legacy_beta_tors"], "1/100")
        self.assertEqual(audit["parameters"]["strict_eta"], "9/5")

    def test_two_point_exact_match_rejects_strict_eta(self):
        matrix = self.payload["damping_compression_bridge_atom_audit"]["two_point_exact_match_matrix"]
        self.assertTrue(matrix["all_pairs_reject_strict_eta_as_legacy_linear_exact_match"])
        for row in matrix["rows"]:
            self.assertEqual(row["forced_eta_for_exact_legacy_linear_match"], "1")
            self.assertFalse(row["strict_eta_matches_forced_eta"])

    def test_effective_beta_varies_and_no_bridge_export(self):
        audit = self.payload["damping_compression_bridge_atom_audit"]
        self.assertFalse(audit["beta_eff_stats"]["constant_across_distances"])
        self.assertGreater(audit["beta_eff_stats"]["max_over_min"], 1.0)
        premise = audit["premise_matrix"]
        self.assertFalse(premise["accepted_as_damping_compression_bridge_atom"])
        self.assertIn("eta_source_law", premise["missing_premises"])
        self.assertIn("target_independent_positive_beta_source", premise["missing_premises"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["full_kernel_bridge_exported"])
        self.assertFalse(flags["role_transfer_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_damping_compression_bridge_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_damping_compression_bridge_atom"])
        self.assertIn("P2849/S1799", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2849/S1799", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2849/S1799", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2849", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
